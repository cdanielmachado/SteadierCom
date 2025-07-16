#!/usr/bin/env python

import argparse
import textwrap
import glob
import os
from math import inf

import pandas as pd
from reframed.community.model import Community
from steadiercom import SteadierCom, SteadierSample
from reframed.solvers.solution import Status
from reframed.io.cache import ModelCache
from reframed import Environment, set_default_solver


def extract_id_from_filepath(filepath):
    filename = os.path.basename(filepath)
    if filename.endswith('.xml'):
        return filename[:-4]
    elif filename.endswith('.xml.gz'):
        return filename[:-7]
    else:
        raise IOError(f'Unrecognized extension in file {filename}. Valid extensions are .xml and .xml.gz')


def build_cache(models):
    ids = [extract_id_from_filepath(model) for model in models]
    return ModelCache(ids, models, load_args={'flavor': 'bigg'})


def load_communities(models, communities):
    if len(models) == 1 and '*' in models[0]:
        models = glob.glob(models[0])
        if not models:
            raise RuntimeError(f'No files found: {models[0]}')
    model_cache = build_cache(models)

    if communities is None:
        return model_cache, {'all': model_cache.get_ids()}, False

    df = pd.read_csv(communities, sep='\t', header=None)
    if len(df.columns) == 2:
        comm_dict = {name: group[1].tolist() for name, group in df.groupby(0)}
        return model_cache, comm_dict, False
    elif len(df.columns) == 3:
        comm_dict = {name: dict(group[[1, 2]].values) for name, group in df.groupby(0)}
        return model_cache, comm_dict, True
    else:
        raise IOError(f'Unexpected number of columns in {communities}')


def load_media_db(filename):
    data = pd.read_csv(filename, sep='\t')
    has_bounds = 'bound' in data.columns
    media_db = {}
    for medium, data_i in data.groupby('medium'):
        if has_bounds:
            media_db[medium] = dict(data_i[['compound', 'bound']].values)
        else:
            media_db[medium] = list(data_i['compound'])
    return media_db, has_bounds


def get_fmt_func(model_format):
    if model_format == 'gapseq':
        return lambda x: f"R_EX_{x}_e0"
    elif model_format == 'bigg':
        return lambda x: f"R_EX_{x}_e"
    else:
        raise ValueError(f"Unsupported model format: {model_format}")


def main_run(models, communities=None, output=None, media=None, mediadb=None, growth=None, sample=None,
             w_e=None, w_r=None, target=None, unlimited=None, model_format='bigg'):

    abstol = 1e-9
    default_growth = 0.1
    fmt_func = get_fmt_func(model_format)

    model_cache, comm_dict, has_abundance = load_communities(models, communities)

    media = media.split(',') if media else [None]
    if any(media) and not mediadb:
        raise RuntimeError('Media database file must be provided.')
    media_db, media_has_bounds = load_media_db(mediadb) if mediadb else ({}, False)

    unlimited_ids = set()
    if unlimited:
        tmp = pd.read_csv(unlimited, header=None)
        unlimited = set(tmp[0])
        suffix = '_e0' if model_format == 'gapseq' else '_e'
        unlimited_ids = {f'M_{x}{suffix}' for x in unlimited}

    if not has_abundance and growth is None:
        growth = default_growth

    results = []

    for comm_id, organisms in comm_dict.items():
        abundance = organisms if has_abundance else None
        org_ids = organisms if not has_abundance else organisms.keys()
        comm_models = [model_cache.get_model(org_id, reset_id=False) for org_id in org_ids]
        community = Community(comm_id, comm_models, copy_models=False)

        # --- Save merged community model for inspection ---
        from reframed.io.sbml import save_model
        output_path = f"{comm_id}_merged_model.xml"
        save_model(community.merged_model, output_path)
        print(f"Saved merged model: {output_path}")

        if model_format == 'gapseq':
            for org_id, org in community.organisms.items():
                org.biomass_reaction = 'R_bio1'

        if target and target not in community.merged_model.reactions:
            raise RuntimeError(f'Invalid target reaction: {target}')

        for medium in media:
            if medium is None:
                medium = 'complete'
                env = Environment.complete(community.merged_model, inplace=False)
            else:
                env = Environment.from_compounds(media_db[medium]).apply(
                    community.merged_model, inplace=False, exclusive=True, warning=False
                )
                if media_has_bounds:
                    for cpd, bound in media_db[medium].items():
                        r_id = f'R_EX_{cpd}_e0' if model_format == 'gapseq' else f'R_EX_{cpd}_e'
                        if r_id in env:
                            env[r_id] = (-bound, inf)

            print(f'Simulating {comm_id} in {medium} medium')

            if unlimited:
                env.update(Environment.from_compounds(unlimited, fmt_func=fmt_func, max_uptake=1000).apply(
                    community.merged_model, inplace=False, exclusive=False, warning=False
                ))

            if sample is None:
                sol = SteadierCom(community, abundance=abundance, growth=growth, allocation=True,
                                  constraints=env, w_e=w_e, w_r=w_r, objective=target)
                if sol.status == Status.OPTIMAL:
                    df = sol.cross_feeding(as_df=True).fillna('environment')
                    df['community'] = comm_id
                    df['medium'] = medium
                    results.append(df)
            else:
                sols = SteadierSample(community, n=sample, abundance=abundance, growth=growth, allocation=True,
                                      constraints=env, w_e=w_e, w_r=w_r, objective=target)
                feasible = [sol.cross_feeding(as_df=True).fillna('environment') for sol in sols if sol.status == Status.OPTIMAL]
                if feasible:
                    df = pd.concat(feasible)
                    df['frequency'] = 1
                    df = df.groupby(['donor', 'receiver', 'compound'], as_index=False).agg({
                        'mass_rate': 'mean',
                        'rate': 'mean',
                        'frequency': lambda x: sum(x) / len(feasible)
                    })
                    df['community'] = comm_id
                    df['medium'] = medium
                    results.append(df)

    if results:
        df_all = pd.concat(results).query(f'rate > {abstol}').sort_values(['community', 'medium', 'mass_rate'], ascending=False)
        if unlimited_ids:
            df_all = df_all[~df_all['compound'].isin(unlimited_ids)]
        output_file = f'{output}.tsv' if output else 'output.tsv'
        df_all.to_csv(output_file, sep='\t', index=False)
        return df_all
    else:
        print('No feasible solutions found.')


def main():
    parser = argparse.ArgumentParser(
        description='Run SteadierCom with BiGG or Gapseq models',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('--format', choices=['gapseq', 'bigg'], default='bigg',
                        help='Specify the model format (e.g., gapseq, bigg)')
    parser.add_argument('models', metavar='MODELS', nargs='+',
                        help="Multiple single-species models (use wildcards like models/*.xml)")
    parser.add_argument('-c', '--communities', metavar='COMMUNITIES.TSV', dest='communities', help="Community definition file")
    parser.add_argument('-o', '--output', dest='output', help="Prefix for output file(s)")
    parser.add_argument('-m', '--media', dest='media', help="Specify a growth medium")
    parser.add_argument('--mediadb', help="Media database (TSV file)")
    parser.add_argument('--growth', type=float, help="Community growth rate (optional)")
    parser.add_argument('--sample', type=int, help="Run sampling analysis for each simulation with N samples")
    parser.add_argument('--we', type=float, default=0.002, help="Weighting factor for enzyme sector")
    parser.add_argument('--wr', type=float, default=0.2, help="Weighting factor for ribosome sector")
    parser.add_argument('--target', help="Target reaction to maximize (optional)")
    parser.add_argument('--solver', help="Select LP solver (options: gurobi, cplex, scip)")
    parser.add_argument('--unlimited', help="Compounds to be considered in excess supply")

    args = parser.parse_args()

    if args.format == 'gapseq':
        print("Using Gapseq format for models")
    else:
        print("Using BiGG format for models")

    if args.solver:
        set_default_solver(args.solver)

    main_run(
        models=args.models,
        communities=args.communities,
        output=args.output,
        media=args.media,
        mediadb=args.mediadb,
        growth=args.growth,
        sample=args.sample,
        w_e=args.we,
        w_r=args.wr,
        target=args.target,
        unlimited=args.unlimited,
        model_format=args.format
    )


if __name__ == '__main__':
    main()


