SPP
===

Codes for analyzing the result of simulation. Contains opengl visualization program and analyzer.

To compile the visual program:
g++ -O3 -lGL -lGLU -lglut -lgsl -lcblas -lopencv_core -lopencv_imgproc -lopencv_calib3d -lopencv_video -lopencv_highgui opengl.cpp -o show
g++ -O3 -lboost_system -lboost_iostreams -lGL -lGLU -lglut -lgsl -lcblas -lopencv_core -lopencv_imgproc -lopencv_calib3d -lopencv_video -lopencv_highgui opengl.cpp -o show

then copy "show" to serial folder 

To compile the analyzer:
g++ -O3 /SPP/analyze/analyze.cpp -lboost_system -lboost_iostreams -lgsl -lcblas -o analyze.out
