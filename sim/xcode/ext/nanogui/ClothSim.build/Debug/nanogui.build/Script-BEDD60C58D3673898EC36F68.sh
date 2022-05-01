#!/bin/sh
set -e
if test "$CONFIGURATION" = "Debug"; then :
  cd /Users/seenum/cs184/Bubble_Rings/sim/xcode/ext/nanogui
  /usr/local/Cellar/cmake/3.22.2/bin/cmake -E cmake_symlink_library /Users/seenum/cs184/Bubble_Rings/sim/xcode/ext/nanogui/Debug/libnanogui.dylib /Users/seenum/cs184/Bubble_Rings/sim/xcode/ext/nanogui/Debug/libnanogui.dylib /Users/seenum/cs184/Bubble_Rings/sim/xcode/ext/nanogui/Debug/libnanogui.dylib
fi
if test "$CONFIGURATION" = "Release"; then :
  cd /Users/seenum/cs184/Bubble_Rings/sim/xcode/ext/nanogui
  /usr/local/Cellar/cmake/3.22.2/bin/cmake -E cmake_symlink_library /Users/seenum/cs184/Bubble_Rings/sim/xcode/ext/nanogui/Release/libnanogui.dylib /Users/seenum/cs184/Bubble_Rings/sim/xcode/ext/nanogui/Release/libnanogui.dylib /Users/seenum/cs184/Bubble_Rings/sim/xcode/ext/nanogui/Release/libnanogui.dylib
fi
if test "$CONFIGURATION" = "MinSizeRel"; then :
  cd /Users/seenum/cs184/Bubble_Rings/sim/xcode/ext/nanogui
  /usr/local/Cellar/cmake/3.22.2/bin/cmake -E cmake_symlink_library /Users/seenum/cs184/Bubble_Rings/sim/xcode/ext/nanogui/MinSizeRel/libnanogui.dylib /Users/seenum/cs184/Bubble_Rings/sim/xcode/ext/nanogui/MinSizeRel/libnanogui.dylib /Users/seenum/cs184/Bubble_Rings/sim/xcode/ext/nanogui/MinSizeRel/libnanogui.dylib
fi
if test "$CONFIGURATION" = "RelWithDebInfo"; then :
  cd /Users/seenum/cs184/Bubble_Rings/sim/xcode/ext/nanogui
  /usr/local/Cellar/cmake/3.22.2/bin/cmake -E cmake_symlink_library /Users/seenum/cs184/Bubble_Rings/sim/xcode/ext/nanogui/RelWithDebInfo/libnanogui.dylib /Users/seenum/cs184/Bubble_Rings/sim/xcode/ext/nanogui/RelWithDebInfo/libnanogui.dylib /Users/seenum/cs184/Bubble_Rings/sim/xcode/ext/nanogui/RelWithDebInfo/libnanogui.dylib
fi

