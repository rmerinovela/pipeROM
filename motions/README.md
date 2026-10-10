# motions/

Floor acceleration histories, described in `motion_sets.yaml`. Format, provided sets and how to add new
ones: [docs/input_files.md](../docs/input_files.md#floor-motions). The floor-motion files (2 GB) are not
tracked by git. They are on Zenodo, [10.5281/zenodo.23283900](https://doi.org/10.5281/zenodo.23283900); `python -m piperom download-motions` downloads them into
`floor_motions/`. The Zenodo record also has the OpenSees model of the 4-storey frame that produced them
and its input ground motions.
