[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![PyPI version](https://badge.fury.io/py/steadiercom.svg)](https://badge.fury.io/py/steadiercom)

![SteadierCom](steadiercom.png)

# SteadierCom
Microbial community simulation using genome-scale metabolic models


## Instalation

```
pip install steadiercom
```

## Usage

Run for a single community without specified media (using complete medium by default)

```
steadiercom species1.xml species2.xml ... speciesN.xml
```

You can also use wildcards: 

```
steadiercom *.xml 
```

You can specify your own growth media in a TSV file. 

```
steadiercom [models] -m M9,LB --mediadb media.tsv
```



You can simulate  multiple communities (specified on a TSV file with the following format **without the header**):

col 1 | col 2
--- | --- 
comm1 | species1
comm1 | species2
comm2 | species1
comm2 | species3


```
steadiercom [models] -c communities.tsv
```

In this situation abundance is variable and growth rate is fixed (default 0.1 h-1). 

You can select a different growth rate:

```
steadiercom [models] [args] --growth 0.2
```

Run for multiple communities with specified relative abundance (**without the header**):


col 1 | col 2 | col 3
--- | --- | ---
comm1 | species1 | 0.5
comm1 | species2 | 0.5
comm2 | species1 | 0.8
comm2 | species3 | 0.2

```
steadiercom [models] -c communities.tsv
```

In this case abundance is fixed and growth rate is maximized (but you can still constrain it with `--growth` if necessary).

To run a sampling analysis and enumerate multiple solutions use:

```
steadiercom [models] [args] --sample 100
```

To change the output filename you can use:

```
steadiercom [models] [args] -o myfile 
```
______ 
## Credits

Daniel Machado, NTNU (2024)

