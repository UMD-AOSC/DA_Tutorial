# Run with input: ETKF, PF, or Hybrid
echo "Usage:"
echo "sh runall_tutorial_4.sh"
echo ""
echo "(Optional)"
echo "sh runall_tutorial_4.sh <method2>"
echo "Hybrid is analyzed and compared to <method2>"

method1=Hybrid
method2=$2
python analysis_init.py $method1
python generate_analysis_3dEns.py
python plot_analysis_vs_nature.py $method1
if [ $# -gt 0 ]; then
  python plot_analysis_vs_analysis.py $method1 $method2
fi
