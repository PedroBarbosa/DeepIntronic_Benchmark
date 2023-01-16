#! /bin/sh
export BUILDAH_FORMAT=docker
docker build . -t intronic_eval
