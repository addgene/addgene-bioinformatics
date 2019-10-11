DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd | sed s%/src/sh%% )"

git clone git@github.com:jfass/apc.git

pushd apc
cp apc.pl $DIR/local-tmp/bin
popd

rm -rf apc
