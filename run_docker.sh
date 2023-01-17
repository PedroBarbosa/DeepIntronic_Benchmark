#! /bin/sh
mkdir -p out
chmod -R 777 out data

docker run -v $PWD/data:/work/data -v $PWD/scripts:/work/scripts -v $PWD/out:/work/out --rm --user=root intronic_eval papermill -k python3 intronic_benchmark.ipynb /work/out/out.ipynb
docker run -v $PWD/data:/work/data -v $PWD/scripts:/work/scripts -v $PWD/out:/work/out --rm --user=root intronic_eval jupyter nbconvert --to html /work/out/out.ipynb
