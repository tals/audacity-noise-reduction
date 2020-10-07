ROOT="$( cd "$(dirname $(dirname "$0"))" >/dev/null 2>&1 ; pwd -P )"
TAG=audiocity-noise-cancellation
pushd $ROOT
docker build -f ci/Dockerfile -t $TAG .
docker run --rm -v $ROOT:/app $TAG "cd app && make && ./bin/noisereduction_driver"