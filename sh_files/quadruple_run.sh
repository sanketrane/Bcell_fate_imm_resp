
while getopts m: flag
do
    case "${flag}" in
        m) modelname=${OPTARG};;
    esac
done

for i in 1 2 3 4
    do
      ./stan_models/ki67_SHM_${modelname} sample num_warmup=300 num_samples=500 data file=datafiles/MZ_data.Rdump \
      output file=save_csv/shm_${modelname}_${i}.csv &
    done
