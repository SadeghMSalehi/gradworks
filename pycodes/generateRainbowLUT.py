#!/usr/bin/python

# from Red (1 0 0) green increases
# when RG (1 1 0) red decreases
# when Green (0 1 0) blue increases
# when GB (0 1 1) green decreases
# when B (0 0 1) red increases
# when RB (1 0 1)

import sys

minIndex = 0
maxIndex = 1900
segment  = (maxIndex - minIndex + 1) / 5.0
reverse = 1

def clamp(x, l, h):
  if (x < l):
    return l
  elif (x > h):
    return h
  return x

def printRGB(idx, rgb):
  print str(idx) + " " + str(idx) + " " + str(int(round(255*rgb[0]))) + " " + str(int(round(255*rgb[1]))) + " " + str(int(round(255*rgb[2]))) + " 255"

def rainbow():
  idx = 0
  if (reverse == 0):
    rgb = [1, 0, 0]
    for y in range(0, 5):
      for x in range(int(round(y*segment)), int(round((y+1)*segment))):
        step = 1 / segment
        if y == 0:
          printRGB(idx, rgb)
          rgb[1] += step
        elif y == 1:
          printRGB(idx, rgb)
          rgb[0] -= step
        elif y == 2:
          printRGB(idx, rgb)
          rgb[2] += step
        elif y == 3:
          printRGB(idx, rgb)
          rgb[1] -= step
        elif y == 4:
          printRGB(idx, rgb)
          rgb[0] += step
        for z in range(0,3):
          rgb[z] = clamp(rgb[z], 0, 1)
        idx += 1
  else:
    rgb = [1, 0, 1]
    for y in range(0, 5):
      for x in range(int(round(y*segment)), int(round((y+1)*segment))):
        step = 1 / segment
        if y == 0:
          printRGB(idx, rgb)
          rgb[0] -= step
        elif y == 1:
          printRGB(idx, rgb)
          rgb[1] += step
        elif y == 2:
          printRGB(idx, rgb)
          rgb[2] -= step
        elif y == 3:
          printRGB(idx, rgb)
          rgb[0] += step
        elif y == 4:
          printRGB(idx, rgb)
          rgb[1] -= step
        for z in range(0,3):
          rgb[z] = clamp(rgb[z], 0, 1)
        idx += 1

rainbow()
