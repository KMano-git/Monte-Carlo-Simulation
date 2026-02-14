# --------------------------------------------------
# make dumy nDdata
# 1st argument:x range, integer
# 2nd argument:y range, integer
# 3rd argument:z range, integer
# --------------------------------------------------
def makeDumy(x, y, z):
    data = []
    # 1D
    if x > 0 and y == 0 and z == 0:
        for i in range(x):
            data.append(0)
    # 2D
    elif x > 0 and y > 0 and z == 0:
        for i in range(y):
            data.append([])
            for j in range(x):
                data[i].append(0)
    # 3D
    elif x > 0 and y > 0 and z > 0:
        for i in range(z):
            data.append([])
            for j in range(y):
                data[i].append([])
                for k in range(x):
                    data[i][j].append(0)

    # return dumy nD array
    return data

if __name__ == "__main__":
    a = makeDumy(6500, 0, 0)
    b = makeDumy(160, 80, 0)
    c = makeDumy(160, 80, 5)
