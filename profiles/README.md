# Execution profiles

`local/` is the conservative workstation default: four cores and a 32 GB
aggregate memory ceiling. `slurm/` enables up to 128 submitted jobs while using
the CPU, memory, disk, and runtime declared by each rule.

Site-specific SLURM accounts, partitions, reservations, and quality-of-service
settings belong in a private overlay or command-line arguments. They are not
hard-coded into the reproducibility archive.
