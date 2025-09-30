        g++ -std=c++17 \
        -Iinclude \
        -I/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/ \
        -I/opt/homebrew/include/ \
        -I/Library/Developer/CommandLineTools/Library/Frameworks/Python3.framework/Versions/3.9/Headers \
        -I/Users/subhodeep/Desktop/CODE/env/lib/python3.10/site-packages/numpy/_core/include \
        src/main.cpp src/RobotDynamics.cpp src/Controller.cpp src/PathPlanningFrenetFrame.cpp src/Path.cpp \
        -L/opt/homebrew/lib \
        -L/Library/Developer/CommandLineTools/Library/Frameworks/Python3.framework/Versions/3.9/lib \
        -lnlopt \
        -lpython3.9 \
        -o main