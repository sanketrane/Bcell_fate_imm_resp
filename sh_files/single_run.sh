
while getopts m: flag
do
    case "${flag}" in
        m) modelname=${OPTARG};;
    esac
done

./stan_models/${modelname} sample num_warmup=500 num_samples=1500 data file=datafiles/Bcell_Imm.Rdump output file=save_csv/sf_${modelname}.csv 
