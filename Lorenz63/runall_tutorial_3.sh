# Run with input: ETKF, PF, or Hybrid
if [ $# -lt 1 ]; then
  echo "Usage:"
  echo "sh runall_tutorial_3.sh <method>"
  echo ""
  echo "(Optional)"
  echo "sh runall_tutorial_3.sh <method1> <method2>"
  echo "<method1> is analyzed and compared to <method2>"
  exit 1
fi
method1=$1
method2=$2
python analysis_init.py $method1
python generate_analysis_3dEns.py
python plot_analysis_vs_nature.py $method1
if [ $# -gt 1 ]; then
  python plot_analysis_vs_analysis.py $method1 $method2
fi
