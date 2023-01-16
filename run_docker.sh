#! /bin/sh
mkdir -p out
chmod -R 777 out

docker run -v $PWD/data:/work/data -v $PWD/scripts:/work/scripts -v $PWD/out:/work/out --rm intronic_eval papermill -k python3 run_benchmark.ipynb /work/out/out.ipynb
docker run -v $PWD/data:/work/data -v $PWD/scripts:/work/scripts -v $PWD/out:/work/out --rm intronic_eval jupyter nbconvert --to html /work/out/out.ipynb
