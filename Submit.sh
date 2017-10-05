Submit.sh
#!/bin/bash
qsub -q \*@@elves -pe mpi-2 2 -l mem=8G -l h_rt=24:00:00 Script2a.sh
qsub -q \*@@elves -pe mpi-spread 2 -l mem=8G -l h_rt=24:00:00 Script2b.sh
qsub -q \*@@elves -pe mpi-4 4 -l mem=8G -l h_rt=12:00:00 Script4a.sh
qsub -q \*@@elves -pe mpi-spread 4 -l mem=8G -l h_rt=12:00:00 Script4b.sh
qsub -q \*@@elves -pe mpi-8 8 -l mem=8G -l h_rt=12:00:00 Script8a.sh
qsub -q \*@@elves -pe mpi-spread 8 -l mem=8G -l h_rt=12:00:00 Script8b.sh
qsub -q \*@@elves -pe mpi-16 16 -l mem=8G -l h_rt=6:00:00 Script16a.sh
qsub -q \*@@elves -pe mpi-spread 16 -l mem=8G -l h_rt=6:00:00 Script16b.sh
qsub -q \*@@elves -pe mpi-8 32 -l mem=8G -l h_rt=6:00:00 Script32a.sh
qsub -q \*@@elves -pe mpi-16 32 -l mem=8G -l h_rt=6:00:00 Script32b.sh
qsub -q \*@@elves -pe mpi-spread 32 -l mem=8G -l h_rt=6:00:00 Script32c.sh
qsub -q \*@@elves -pe mpi-8 64 -l mem=8G -l h_rt=6:00:00 Script64a.sh
qsub -q \*@@elves -pe mpi-16 64 -l mem=8G -l h_rt=6:00:00 Script64b.sh
qsub -q \*@@elves -pe mpi-spread 64 -l mem=8G -l h_rt=6:00:00 Script64c.sh
qsub -q \*@@elves -pe mpi-2 16 -l mem=8G -l h_rt=6:00:00 Script16c.sh
qsub -q \*@@elves -pe mpi-4 16 -l mem=8G -l h_rt=6:00:00 Script16d.sh
qsub -q \*@@elves -pe mpi-8 16 -l mem=8G -l h_rt=6:00:00 Script16e.sh
