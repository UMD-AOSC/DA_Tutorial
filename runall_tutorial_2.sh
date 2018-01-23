# Run with input: nudging, OI, or 3DVar
echo "Usage:"
echo "sh runall_tutorial_2.sh <method>"
echo "<method1> is analyzed and compared to <method2>"
method1=$1
method2=$2
python analysis_init.py $method1
python generate_analysis_3dDet.py
python plot_analysis_vs_nature.py $method1
python plot_analysis_vs_analysis.py $method1 $method2
