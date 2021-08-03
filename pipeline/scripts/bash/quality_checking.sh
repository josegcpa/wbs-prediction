NET_PATH=scripts/python/attention_classifier.py
#NET_PATH=/nfs/research1/gerstung/josegcpa/projects/01IMAGE/quality_net/attention_classifier.py
CKPT_PATH=/nfs/research1/gerstung/josegcpa/projects/01IMAGE/quality_net/achtung_no_attention.ckpt
SLIDE_PATH=$1
PLACEHOLDER_STR=n
PLACEHOLDER_FLOAT=0.01
PLACEHOLDER_INT=1

python3 $NET_PATH\
  $CKPT_PATH\
  $SLIDE_PATH\
  0.0001\
  100000\
  16\
  predict\
  N\
  False\
  2

#python3 $NET_PATH\
#  $CKPT_PATH\
#  $SLIDE_PATH\
#  128\
#  False\
#  2
