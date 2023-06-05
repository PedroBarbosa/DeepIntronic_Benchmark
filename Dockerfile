FROM docker.io/mambaorg/micromamba:0.27.0
LABEL Name=DeepIntronicBenchmark Version=1.0

COPY conda_env.yml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes
WORKDIR /work
COPY intronic_benchmark.ipynb /work/intronic_benchmark.ipynb
VOLUME [ "/work/data" ]
VOLUME [ "/work/scripts" ]
VOLUME [ "/work/out" ]
