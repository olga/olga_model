import sys

if (len(sys.argv) == 2):
    filename = sys.argv[1]
else:
    raise RuntimeError("Illegal number of arguments")

fileread = open(filename, "r")
lines = fileread.readlines()
fileread.close()

# Find the line number where the copyright info starts.
nline = 0
for n in lines:
  # Store the commenting style, to make sure C++, Python and Cmake work.
  npos = n.find('Copyright')
  if (npos != -1):
    left_of_copyright = n[0:npos]
    break
  nline += 1

# Delete the three lines of this header.
del(lines[nline:nline+1])

newlines = ['{0}Copyright (c) 2013-2016 Bart van Stratum\n'.format(left_of_copyright),
            '{0}Copyright (c) 2015-2016 Roel Baardman\n'.format(left_of_copyright)]

# Insert the new header.
lines[nline:nline] = newlines[:]

# Save the output.
filewrite = open(filename, "w")
filewrite.writelines(lines)
filewrite.close()
