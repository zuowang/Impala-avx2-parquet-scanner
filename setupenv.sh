#!/usr/bin/env bash

export JAVA_HOME=/usr/lib/jvm/java-7-oracle
export ANT_OPTS="-Dhttp.proxyHost=child-prc.intel.com -Dhttp.proxyPort=913"
export GIT_PROXY_COMMAND=/usr/local/bin/socks-gateway

export IMPALA_HOME=/mnt/disk1/dev/Impala
export PIC_LIB_PATH=${IMPALA_HOME}
#export HADOOP_LZO=/mnt/disk1/dev/Impala/../hadoop-lzo
#export IMPALA_LZO=/mnt/disk1/dev/Impala/../impala-lzo

#. ${IMPALA_HOME}/impala-config.sh
export PRE_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
. ${IMPALA_HOME}/bin/impala-config.sh
. ${IMPALA_HOME}/bin/set-classpath.sh
. ${IMPALA_HOME}/bin/set-pythonpath.sh
export LD_LIBRARY_PATH=$PRE_LD_LIBRARY_PATH
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
