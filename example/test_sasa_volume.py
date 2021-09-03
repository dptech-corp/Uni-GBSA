# from ..lib.pdbremix import 
import ..pdbremix as pdbremix
def main():
    pdbfile = './2fvy.pdb'
    vol = pdbremix.cal_volume(pdbfile)
    print(vol)
    sas = pdbremix.cal_sasa(pdbfile)
    print(sas)

if __name__ == "__main__":
    main()