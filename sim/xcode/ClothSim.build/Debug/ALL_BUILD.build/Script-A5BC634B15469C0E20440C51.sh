#!/bin/sh
set -e
if test "$CONFIGURATION" = "Debug"; then :
  cd /Users/seenum/cs184/Bubble_Rings/sim/xcode
  echo Build\ all\ projects
fi
if test "$CONFIGURATION" = "Release"; then :
  cd /Users/seenum/cs184/Bubble_Rings/sim/xcode
  echo Build\ all\ projects
fi
if test "$CONFIGURATION" = "MinSizeRel"; then :
  cd /Users/seenum/cs184/Bubble_Rings/sim/xcode
  echo Build\ all\ projects
fi
if test "$CONFIGURATION" = "RelWithDebInfo"; then :
  cd /Users/seenum/cs184/Bubble_Rings/sim/xcode
  echo Build\ all\ projects
fi

