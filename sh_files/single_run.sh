
while getopts m: flag
do
    case "${flag}" in
        m) modelname=${OPTARG};;
    esac
done

stan_models/MAP_${modelname} sample num_samples=500 num_warmup=300 data file=data/Treg_data_S33.Rdump output file=save_csv/${modelname}_c1.csv