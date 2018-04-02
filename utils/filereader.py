def readLinesTo(filename, startObject, lineFunc):
    f = open(filename, "r")
    for line in f:
        startObject = lineFunc(line, startObject)
    return startObject
