
while getopts m: flag
do
    case "${flag}" in
        m) modelname=${OPTARG};;
    esac
done

./stan_models/New_modelsCAR/${modelname} sample num_warmup=500 num_samples=500 data file=datafiles/Bcell_Imm.Rdump output file=save_csv/${modelname}_${i}.csv
