#!/usr/bin/env python

import argparse
import textwrap
from reframed.community.model import Community
from .steadiercom import SteadierCom, SteadierSample
from reframed.solvers.solution import Status

import glob
import pandas as pd
from reframed.io.cache import ModelCache
import os
from reframed import Environment, set_default_solver


def extract_id_from_filepath(filepath):
    filename = os.path.basename(filepath)

    if filename.endswith('.xml'):
        organism_id = filename[:-4]
    elif filename.endswith('.xml.gz'):
        organism_id = filename[:-7]
    else:
        raise IOError(f'Unrecognized extension in file {filename}. Valid extensions are .xml and .xml.gz')

    return organism_id


def build_cache(models):
    ids = [extract_id_from_filepath(model) for model in models]

    return ModelCache(ids, models, load_args={'flavor': 'bigg'})


def load_communities(models, communities):
    if len(models) == 1 and '*' in models[0]:
        pattern = models[0]
        models = glob.glob(pattern)
        if len(models) == 0:
            raise RuntimeError(f'No files found: {pattern}')
        
    model_cache = build_cache(models)

    has_abundance = False

    if communities is None:
        comm_dict = {'all': model_cache.get_ids()}
    else:
        df = pd.read_csv(communities, sep='\t', header=None)

        if len(df.columns) == 2:
            comm_dict = {name: group[1].tolist() for name, group in df.groupby(0)}
        elif len(df.columns) == 3:
            comm_dict = {name: dict(group[[1,2]].values) for name, group in df.groupby(0)}
            has_abundance = True
        else:
            raise IOError(f'Unexpected number of columns in {communities}')

    return model_cache, comm_dict, has_abundance


def load_media_db(filename, sep='\t', medium_col='medium', compound_col='compound'):
    """ Load media library file. """

    data = pd.read_csv(filename, sep=sep)
    media_db = data[[medium_col, compound_col]].groupby(medium_col).agg(lambda x: list(x))

    return media_db[compound_col].to_dict()


def main_run(models, communities=None, output=None, media=None, mediadb=None, growth=None, sample=None, w_e=None, w_r=None, target=None):

    abstol = 1e-9
    default_growth = 0.1

    model_cache, comm_dict, has_abundance = load_communities(models, communities)

    if media is None:
        media = [None]
    else:
        media = media.split(',')

        if mediadb is None:
            raise RuntimeError('Media database file must be provided.')
        else:   
            media_db = load_media_db(mediadb)

    results = []
    
    if not has_abundance and growth is None:
        growth = default_growth

    for comm_id, organisms in comm_dict.items():

        if has_abundance:
            abundance = organisms
            organisms = organisms.keys()
        else:
            abundance = None

        comm_models = [model_cache.get_model(org_id, reset_id=True) for org_id in organisms]
        community = Community(comm_id, comm_models, copy_models=False)

        if target is not None and target not in community.merged_model.reactions:
            raise RuntimeError(f'Invalid target reaction: {target}')

        for medium in media:

            if medium is None:
                print(f'simulating {comm_id} in complete medium')
                env = Environment.complete(community.merged_model, inplace=False)
            else:
                print(f'simulating {comm_id} in {medium}')
                env = Environment.from_compounds(media_db[medium]).apply(community.merged_model, inplace=False, exclusive=True, warning=False)

            if sample is None:
                sol = SteadierCom(community, abundance=abundance, growth=growth, allocation=True, constraints=env, w_e=w_e, w_r=w_r, objective=target)
                if sol.status == Status.OPTIMAL:
                    df = sol.cross_feeding(as_df=True).fillna('environment')
                    df['community'] = comm_id
                    df['medium'] = medium
                    results.append(df)
            else:
                sols = SteadierSample(community, n=sample, abundance=abundance, growth=growth, allocation=True, constraints=env, w_e=w_e, w_r=w_r, objective=target)
                feasible = [sol.cross_feeding(as_df=True).fillna('environment') for sol in sols if sol.status == Status.OPTIMAL]
                if len(feasible) > 0:
                    df = pd.concat(feasible)
                    df['frequency'] = 1
                    df = df.groupby(['donor', 'receiver', 'compound'], as_index=False).agg(
                        {'mass_rate': lambda x: x.mean(), 'rate': lambda x: x.mean(), 'frequency': lambda x: sum(x)/len(feasible)})
                    df['community'] = comm_id
                    df['medium'] = medium
                    results.append(df)

    if not output:
        output_file = 'output.tsv'
    else:
        output_file = f'{output}.tsv'

    if len(results) > 0:
        df_all = pd.concat(results).query(f'rate > {abstol}').sort_values('mass_rate', ascending=False)
        df_all.to_csv(output_file, sep='\t', index=False)
        return df_all
    else:
        print('No feasible solutions found.')


def main():

    parser = argparse.ArgumentParser(description="Simulate microbial communities with SteadierCom",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('models', metavar='MODELS', nargs='+',
                        help=textwrap.dedent(
        """
        Multiple single-species models (one or more files).
        
        You can use wild-cards, for example: models/*.xml, and optionally protect with quotes to avoid automatic bash
        expansion (this will be faster for long lists): "models/*.xml". 
        """
        ))

    parser.add_argument('-c', '--communities', metavar='COMMUNITIES.TSV', dest='communities',
                        help=textwrap.dedent(
        """
        Run SteadierCom for multiple (sub)communities.
        The communities must be specified in a tab-separated file with community and organism identifiers.
        The organism identifiers should match the file names in the SBML files (without extension).
        
        Example:
            community1\torganism1
            community1\torganism2
            community2\torganism1
            community2\torganism3

        You can optionally specify the relative abundances for each community member:
        
        Example:
            community1\torganism1\t0.5
            community1\torganism2\t0.5
            community2\torganism1\t0.9
            community2\torganism3\t0.1

        """
    ))

    parser.add_argument('-o', '--output', dest='output', help="Prefix for output file(s)")
    parser.add_argument('-m', '--media', dest='media', help="Specify a growth medium")
    parser.add_argument('--mediadb', help="Media database file")
    parser.add_argument('--growth', type=float, help="Community growth rate (optional)")
    parser.add_argument('--sample', type=int, help="Run sampling analysis for each simulation with N samples")
    parser.add_argument('--we', type=float, default=0.002, help="Weighting factor for enzyme sector (default: 0.002 gDW.h/mmol)")
    parser.add_argument('--wr', type=float, default=0.2, help="Weighting factor for ribosome sector (default: 0.2 h)")
    parser.add_argument('--target', help="Target reaction to maximize (optional)")
    parser.add_argument('--solver', help="Select LP solver (options: gurobi, cplex, scip)")

    args = parser.parse_args()

    if args.solver is not None:
        set_default_solver(args.solver)

    main_run(
        models=args.models,
        communities=args.communities,
        output=args.output,
        media=args.media,
        mediadb=args.mediadb,
        growth=args.growth,
        sample=args.sample,
        w_e = args.we,
        w_r = args.wr,
        target = args.target,
    )


if __name__ == '__main__':
    main()