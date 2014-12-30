SPP
===

Serial code.

To compile the simulator:
g++ -lgsl -lcblas -O3 main.cpp

To compile the visual program:
g++ -O3 -lGL -lGLU -lglut -lgsl -lcblas -lopencv_core -lopencv_imgproc -lopencv_calib3d -lopencv_video -lopencv_highgui opengl.cpp -o show
g++ -O3 -lboost_system -lboost_iostreams -lGL -lGLU -lglut -lgsl -lcblas -lopencv_core -lopencv_imgproc -lopencv_calib3d -lopencv_video -lopencv_highgui opengl.cpp -o show

To compile the analyzer:
g++ -O3 analyze.cpp -lboost_system -lboost_iostreams -lGL -lGLU -lglut -lgsl -lcblas -o analyze.out
