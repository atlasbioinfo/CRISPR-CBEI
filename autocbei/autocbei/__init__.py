
def mainCBEI():
    from sys import version_info
    if version_info.major != 3:
        print("Warning!")
        print('\tYour python version is '+str(version_info.major)+'.'+str(version_info.minor))
        print('\t"autocbei" relies on python3, please use python3 (Recommended python>=3.6.0)!')
        print('autocbei exit...')
        return
    from autocbei.cbei.autoCBEI import mainAutoCBEI

    mainAutoCBEI()


if __name__ == "__main__":
    mainCBEI()