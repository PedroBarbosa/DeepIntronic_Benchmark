#! /bin/sh
export BUILDAH_FORMAT=docker
docker build . --no-cache -t intronic_eval
