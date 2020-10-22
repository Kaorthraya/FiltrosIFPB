import numpy as np

def printTF(num, den):
    div = '-----------'
    print('\n')
    num1 = ""
    
    for i in range(0, len(num)):
        if(num[i] != 0):
            if(i != len(num)-1):
                num1 = num1 + "(" + "{:.3e}".format(num[i]) + ")" + '*s^' + str(len(num)-i-1) + " + "
            else:
                num1 = num1 + "(" + "{:.3e}".format(num[-1]) + ")"
    if(num1[-2] == "+"):
        num1 = num1[:len(num1)-2]

    for i in range(1, len(den)):
        div = div + '------------------'

    print(num1.center(len(div)))
    print(div)

    for i in range(0, len(den)):   
        if(i != len(den)-1):
            print("(" + "{:.3e}".format(den[i]) + ")" + '*s^' + str(len(den)-i-1) + " + ", end = "")
        else:
            print("(" + "{:.3e}".format(den[-1]) + ")")
    print('\n')
