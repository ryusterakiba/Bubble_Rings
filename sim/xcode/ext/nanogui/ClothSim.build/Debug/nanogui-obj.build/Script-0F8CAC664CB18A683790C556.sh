#!/bin/sh
set -e
if test "$CONFIGURATION" = "Debug"; then :
  cd /Users/seenum/cs184/Bubble_Rings/sim/xcode/ext/nanogui
  /usr/local/Cellar/cmake/3.22.2/bin/cmake -DOUTPUT_C=nanogui_resources.cpp -DOUTPUT_H=nanogui_resources.h -DINPUT_FILES=/Users/seenum/cs184/Bubble_Rings/sim/ext/nanogui/resources/Roboto-Bold.ttf,/Users/seenum/cs184/Bubble_Rings/sim/ext/nanogui/resources/Roboto-Regular.ttf,/Users/seenum/cs184/Bubble_Rings/sim/ext/nanogui/resources/entypo.ttf -P /Users/seenum/cs184/Bubble_Rings/sim/ext/nanogui/resources/bin2c.cmake
fi
if test "$CONFIGURATION" = "Release"; then :
  cd /Users/seenum/cs184/Bubble_Rings/sim/xcode/ext/nanogui
  /usr/local/Cellar/cmake/3.22.2/bin/cmake -DOUTPUT_C=nanogui_resources.cpp -DOUTPUT_H=nanogui_resources.h -DINPUT_FILES=/Users/seenum/cs184/Bubble_Rings/sim/ext/nanogui/resources/Roboto-Bold.ttf,/Users/seenum/cs184/Bubble_Rings/sim/ext/nanogui/resources/Roboto-Regular.ttf,/Users/seenum/cs184/Bubble_Rings/sim/ext/nanogui/resources/entypo.ttf -P /Users/seenum/cs184/Bubble_Rings/sim/ext/nanogui/resources/bin2c.cmake
fi
if test "$CONFIGURATION" = "MinSizeRel"; then :
  cd /Users/seenum/cs184/Bubble_Rings/sim/xcode/ext/nanogui
  /usr/local/Cellar/cmake/3.22.2/bin/cmake -DOUTPUT_C=nanogui_resources.cpp -DOUTPUT_H=nanogui_resources.h -DINPUT_FILES=/Users/seenum/cs184/Bubble_Rings/sim/ext/nanogui/resources/Roboto-Bold.ttf,/Users/seenum/cs184/Bubble_Rings/sim/ext/nanogui/resources/Roboto-Regular.ttf,/Users/seenum/cs184/Bubble_Rings/sim/ext/nanogui/resources/entypo.ttf -P /Users/seenum/cs184/Bubble_Rings/sim/ext/nanogui/resources/bin2c.cmake
fi
if test "$CONFIGURATION" = "RelWithDebInfo"; then :
  cd /Users/seenum/cs184/Bubble_Rings/sim/xcode/ext/nanogui
  /usr/local/Cellar/cmake/3.22.2/bin/cmake -DOUTPUT_C=nanogui_resources.cpp -DOUTPUT_H=nanogui_resources.h -DINPUT_FILES=/Users/seenum/cs184/Bubble_Rings/sim/ext/nanogui/resources/Roboto-Bold.ttf,/Users/seenum/cs184/Bubble_Rings/sim/ext/nanogui/resources/Roboto-Regular.ttf,/Users/seenum/cs184/Bubble_Rings/sim/ext/nanogui/resources/entypo.ttf -P /Users/seenum/cs184/Bubble_Rings/sim/ext/nanogui/resources/bin2c.cmake
fi

