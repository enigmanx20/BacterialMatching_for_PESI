# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import scipy as sp
from StringIO import StringIO

class PESI_single_data(object):
    """
    data loading and vextorization
    m/z compression
    centroid
    peak detection
    intensity normalization
    scaling
    """
    
    def __init__(self, data_name):
        self.name =  data_name
        
    def __str__(self):
        return self.name
    
    def load(self, path):
        """
        load text file in the path and return a dataframe structured by 
    
        Parameters:
            datatype
            1. LabSolutions ASCII format
            'm/z', 'Absolute Intensity', and 'Relative Intensity'
            filling NaN with 0
       
    
        """
                                
        f = open(path,'r')  #fileを開く
        raw_data = f.read() #データを読み出す
        f.close()
    
        #LabSolutions ASCII format
        index_start = raw_data.find("Profile Data")
        raw_data = raw_data[index_start:]
        index_next = raw_data.find('\r\n')
        raw_data = raw_data[index_next+2:]
        index_next = raw_data.find('\r\n')
        raw_data = raw_data[index_next+2:]
        data = pd.read_csv(StringIO(raw_data), sep='\t')
        
        data = data.fillna(value=0)
        data["m/z"] = np.ceil(data["m/z"] * 10.0) * 0.1  #now m/z starts from 10.1 
            
        self.data = data
        self.path = path
              
    def load_peaks(self, path):
        f = open(path,'rU')  #fileを開く
        raw_data = f.read() #データを読み出す
        f.close()
        
        #LabSolutions ASCII format
        index_start = raw_data.find("Absolute")
        
        index_end = raw_data.find("Profile Data")
        raw_data = raw_data[index_start-4:index_end]
        
        data = pd.read_csv(StringIO(raw_data), sep='\t')
            
        data = data.fillna(value=0)
        data["m/z"] = np.round(data["m/z"] * 10.0) * 0.1  #now m/z starts from 10.1 
                
        self.peak = data
        self.path = path
    
    def load_unitPeaks(self, path):
        f = open(path,'rU')  #fileを開く
        raw_data = f.read() #データを読み出す
        f.close()
        
        #LabSolutions ASCII format
        index_start = raw_data.find("Absolute")
        
        index_end = raw_data.find("Profile Data")
        raw_data = raw_data[index_start-4:index_end]
        
        data = pd.read_csv(StringIO(raw_data), sep='\t')
            
        data = data.fillna(value=0)
        data["m/z"] = np.floor(data["m/z"])  #now m/z starts from 10.1 
        
                
        self.unitPeak = data
        self.path = path
    
    
    def centroid(self):
        self.data["centroid"] = zeros(len(self.data))
        for i  in range(len(self.peak["m/z"])):
                truth_vec = self.data["m/z"] == ones(len(self.data)) * self.peak.iloc[i]["m/z"]
   
                self.data["centroid"][truth_vec] = self.peak.iloc[i]["Relative Intensity"]
            
    
    def compress(self, datatype=1,ctype=1, bin=10):
        """
        要修正
        compress data using "Relative Intensity" 
        Parameters:
        ctype
        1. taking local mean
        
        2. taking local median
        
        3. taking local maxima
        
        4. taking local minima
        
        bin:
        number of bins compressed together
        10 bins  equivelent to each m/z 1.0
        
        """
        start_index = 4  #start from m/z 10.4
        end_index = 19894 #end at m/z 1999.5
        new_intensity = np.zeros(1989)
        new_mz = np.linspace(11,1999,num = 1989)
        label = ["m/z", "New Intensity"]
        
        if ctype == 1: #taking local mean   
            for i in range(1989):
                index_start = 4 + i*10
                new_intensity[i] = np.mean(self.data[index_start:index_start+10]["Relative Intensity"])
            data = pd.DataFrame(np.array([new_mz, new_intensity]).T, columns=["m/z", "New Intensity"])
            self.cdata1 = data
            
        if ctype == 2:#taking local median
            for i in range(1989):
                index_start = 4 + i*10
                new_intensity[i] = np.median(self.data[index_start:index_start+10]["Relative Intensity"])
            data = pd.DataFrame(np.array([new_mz, new_intensity]).T,columns=["m/z", "New Intensity"])
            self.cdata2 = data    
            
        if ctype == 3:#taking local maxima
            for i in range(1989):
                index_start = 4 + i*10
                new_intensity[i] = np.amax(self.data[index_start:index_start+10]["Relative Intensity"])
            data = pd.DataFrame(np.array([new_mz, new_intensity]).T,columns=["m/z", "New Intensity"])
            self.cdata3 = data
            
        if ctype == 4: #taking local minima
            for i in range(1989):
                index_start = 4 + i*10
                new_intensity[i] = np.amin(self.data[index_start:index_start+10]["Relative Intensity"])
            data = pd.DataFrame(np.array([new_mz, new_intensity]).T,columns=["m/z", "New Intensity"])
            self.cdata4 = data
        
        self.cdata = data

"""
Data setting
"""       
dataListP = []
dataListN = []
unknownListP = []
unknownListN = []
bloodP = []
bloodN = []

S_aurP = PESI_single_data("S_aureus_posi"); dataListP.append(S_aurP)
S_aurN = PESI_single_data("S_aureus_nega"); dataListN.append(S_aurN)
S_epiP = PESI_single_data("S_epidermidis_posi"); dataListP.append(S_epiP)
S_epiN = PESI_single_data("S_epidermidis_nega"); dataListN.append(S_epiN)
S_agaP= PESI_single_data("S_agalactiae_posi");dataListP.append(S_agaP)
S_agaN = PESI_single_data("S_agalactiae_nega");dataListN.append(S_agaN)
S_virP = PESI_single_data("S_viridans_posi");dataListP.append(S_virP)
S_virN = PESI_single_data("S_viridans_nega");dataListN.append(S_virN)
#S_faeP = PESI_single_data("S_faecalis_posi"); dataListP.append(S_faeP)
#S_faeN = PESI_single_data("S_faecalis_nega"); dataListN.append(S_faeN)

E_coliP =  PESI_single_data("E_coli_posi"); dataListP.append(E_coliP)
E_coliN =  PESI_single_data("E_coli_nega"); dataListN.append(E_coliN)
P_aerP = PESI_single_data("P_auruginosa_posi"); dataListP.append(P_aerP)
P_aerN = PESI_single_data("P_auruginosa_nega"); dataListN.append(P_aerN)
K_pneP =  PESI_single_data("K_pneumoniae_posi"); dataListP.append(K_pneP)
K_pneN =  PESI_single_data("K_pneumoniae_nega"); dataListN.append(K_pneN)
P_vulP = PESI_single_data("P_vulgaris_posi"); dataListP.append(P_vulP)
P_vulN = PESI_single_data("P_vulgaris_nega"); dataListN.append(P_vulN)
C_albP = PESI_single_data("C_albicans_posi"); dataListP.append(C_albP)
C_albN = PESI_single_data("C_albicans_nega"); dataListN.append(C_albN)



Unknown1P = PESI_single_data("Unknown1P"); unknownListP.append(Unknown1P)
Unknown1N = PESI_single_data("Unknown1N"); unknownListN.append(Unknown1N)
Unknown2P = PESI_single_data("Unknown2P"); unknownListP.append(Unknown2P)
Unknown2N = PESI_single_data("Unknown2N"); unknownListN.append(Unknown2N)
Unknown3P = PESI_single_data("Unknown3P"); unknownListP.append(Unknown3P)
Unknown3N = PESI_single_data("Unknown3N"); unknownListN.append(Unknown3N)
Unknown4P = PESI_single_data("Unknown4P"); unknownListP.append(Unknown4P)
Unknown4N = PESI_single_data("Unknown4N"); unknownListN.append(Unknown4N)

Blood_5perP = PESI_single_data("Blood_5percent_posi"); bloodP.append(Blood_5perP)
Blood_5perN = PESI_single_data("Blood_5percent_nega"); bloodN.append(Blood_5perN)
Blood_S_aurP = PESI_single_data("Blood_S_aureus_posi"); bloodP.append(Blood_S_aurP)
Blood_S_aurN = PESI_single_data("Blood_S_aureus_nega"); bloodN.append(Blood_S_aurN)
Blood_E_coliP = PESI_single_data("Blood_E_coli_posi"); bloodP.append(Blood_E_coliP)
Blood_E_coliN = PESI_single_data("Blood_E_coli_nega"); bloodN.append(Blood_E_coliN)
Blood_P_vulP = PESI_single_data("Blood_P_vulgaris_posi"); bloodP.append(Blood_P_vulP)
Blood_P_vulN = PESI_single_data("Blood_P_vulgaris_nega"); bloodN.append(Blood_P_vulN)


#load data
S_aurP.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140530/S_aureus_posi.txt')
S_aurN.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140530/P_aer_nega.txt')
S_epiP.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140603/S_epi_posi.txt')
S_epiN.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140603/S_epi_nega.txt')
S_agaP.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140706/140706_S_aga_posi.txt')
S_agaN.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140706/140706_S_aga_nega.txt')
S_virP.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140706/140706_S_vir_posi.txt')
S_virN.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140706/140706_S_vir_nega.txt')
#S_faeP.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140612/S_faecalis_posi.txt')
#S_faeN.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140612/S_faecalis_negat.txt')

E_coliP.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140530/E_coli_posi.txt')
E_coliN.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140530/E_coli_nega.txt')
P_aerP.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140530/P_aer_posi.txt')
P_aerN.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140530/P_aer_nega.txt')
K_pneP.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140603/K_pneu_posi.txt')
K_pneN.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140603/K_pneu_nega.txt')
P_vulP.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140612/P_vul_posi.txt')
P_vulN.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140612/P_vul_nega.txt')
C_albP.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/C_alb_posi.txt')
C_albN.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/C_alb_nega.txt')
            

Unknown1P.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Yoshimura_white_posi.txt')
Unknown1N.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Yoshimura_white_nega.txt')
Unknown2P.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Iizuka_white_posi.txt')
Unknown2N.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Iizuka_white_nega.txt')
Unknown3P.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Kawai_white_posi.txt')
Unknown3N.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Kawai_white_negat.txt')
Unknown4P.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Kawai_yellow_posi.txt')
Unknown4N.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Kawai_yellow_nega.txt')

Blood_5perP.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140612/mouse_blood_5per_posi.txt')
Blood_5perN.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140612/mouse_blood_5per_nega.txt')
Blood_S_aurP.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/mouse_blood_S_aur_posi.txt')
Blood_S_aurN.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/mouse_blood_S_aur_nega.txt')
Blood_E_coliP.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/mouse_blood_E_coli_posi.txt')
Blood_E_coliN.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/mouse_blood_E_coli_nega.txt')
Blood_P_vulP.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/mouse_blood_P_vul_posi.txt')
Blood_P_vulN.load('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/mouse_blood_P_vul_nega.txt')

#load peaks
S_aurP.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140530/S_aureus_posi.txt')
S_aurN.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140530/P_aer_nega.txt')
S_epiP.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140603/S_epi_posi.txt')
S_epiN.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140603/S_epi_nega.txt')
S_agaP.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140706/140706_S_aga_posi.txt')
S_agaN.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140706/140706_S_aga_nega.txt')
S_virP.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140706/140706_S_vir_posi.txt')
S_virN.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140706/140706_S_vir_nega.txt')
#S_faeP.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140612/S_faecalis_posi.txt')
#S_faeN.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140612/S_faecalis_negat.txt')

E_coliP.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140530/E_coli_posi.txt')
E_coliN.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140530/E_coli_nega.txt')
P_aerP.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140530/P_aer_posi.txt')
P_aerN.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140530/P_aer_nega.txt')
K_pneP.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140603/K_pneu_posi.txt')
K_pneN.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140603/K_pneu_nega.txt')
P_vulP.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140612/P_vul_posi.txt')
P_vulN.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140612/P_vul_nega.txt')
C_albP.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/C_alb_posi.txt')
C_albN.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/C_alb_nega.txt')

Unknown1P.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Yoshimura_white_posi.txt')
Unknown1N.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Yoshimura_white_nega.txt')
Unknown2P.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Iizuka_white_posi.txt')
Unknown2N.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Iizuka_white_nega.txt')
Unknown3P.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Kawai_white_posi.txt')
Unknown3N.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Kawai_white_negat.txt')
Unknown4P.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Kawai_yellow_posi.txt')
Unknown4N.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Kawai_yellow_nega.txt')

Blood_5perP.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140612/mouse_blood_5per_posi.txt')
Blood_5perN.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140612/mouse_blood_5per_nega.txt')
Blood_S_aurP.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/mouse_blood_S_aur_posi.txt')
Blood_S_aurN.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/mouse_blood_S_aur_nega.txt')
Blood_E_coliP.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/mouse_blood_E_coli_posi.txt')
Blood_E_coliN.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/mouse_blood_E_coli_nega.txt')
Blood_P_vulP.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/mouse_blood_P_vul_posi.txt')
Blood_P_vulN.load_peaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/mouse_blood_P_vul_nega.txt')

#load unit peaks
S_aurP.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140530/S_aureus_posi.txt')
S_aurN.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140530/P_aer_nega.txt')
S_epiP.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140603/S_epi_posi.txt')
S_epiN.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140603/S_epi_nega.txt')
S_agaP.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140706/140706_S_aga_posi.txt')
S_agaN.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140706/140706_S_aga_nega.txt')
S_virP.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140706/140706_S_vir_posi.txt')
S_virN.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140706/140706_S_vir_nega.txt')
#S_faeP.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140612/S_faecalis_posi.txt')
#S_faeN.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140612/S_faecalis_negat.txt')

E_coliP.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140530/E_coli_posi.txt')
E_coliN.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140530/E_coli_nega.txt')
P_aerP.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140530/P_aer_posi.txt')
P_aerN.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140530/P_aer_nega.txt')
K_pneP.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140603/K_pneu_posi.txt')
K_pneN.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140603/K_pneu_nega.txt')
P_vulP.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140612/P_vul_posi.txt')
P_vulN.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140612/P_vul_nega.txt')
C_albP.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/C_alb_posi.txt')
C_albN.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/C_alb_nega.txt')

Unknown1P.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Yoshimura_white_posi.txt')
Unknown1N.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Yoshimura_white_nega.txt')
Unknown2P.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Iizuka_white_posi.txt')
Unknown2N.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Iizuka_white_nega.txt')
Unknown3P.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Kawai_white_posi.txt')
Unknown3N.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Kawai_white_negat.txt')
Unknown4P.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Kawai_yellow_posi.txt')
Unknown4N.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140604/Kawai_yellow_nega.txt')

Blood_5perP.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140612/mouse_blood_5per_posi.txt')
Blood_5perN.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140612/mouse_blood_5per_nega.txt')
Blood_S_aurP.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/mouse_blood_S_aur_posi.txt')
Blood_S_aurN.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/mouse_blood_S_aur_nega.txt')
Blood_E_coliP.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/mouse_blood_E_coli_posi.txt')
Blood_E_coliN.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/mouse_blood_E_coli_nega.txt')
Blood_P_vulP.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/mouse_blood_P_vul_posi.txt')
Blood_P_vulN.load_unitPeaks('/Users/cellbiology/Desktop/PESI_DATA/PESI_MS_sepsis/140613/mouse_blood_P_vul_nega.txt')

#centroid
S_aurP.centroid()
S_aurN.centroid()
S_epiP.centroid()
S_epiN.centroid()
S_agaP.centroid()
S_agaN.centroid()
S_virP.centroid()
S_virN.centroid()
#S_faeP.centroid()
#S_faeN.centroid()

E_coliP.centroid()
E_coliN.centroid()
P_aerP.centroid()
P_aerN.centroid()

Unknown1P.centroid()
Unknown1N.centroid()
Unknown2P.centroid()
Unknown2N.centroid()
Unknown3P.centroid()
Unknown3N.centroid()
Unknown4P.centroid()
Unknown4N.centroid()

Blood_5perP.centroid()
Blood_5perN.centroid()
Blood_S_aurP.centroid()
Blood_S_aurN.centroid()
Blood_E_coliP.centroid()
Blood_E_coliN.centroid()
Blood_P_vulP.centroid()
Blood_P_vulN.centroid()

#define how many pekas to use
numPeaks = 40

"""""
Unknown data matching
"""""
def UnknownMatching(index):
    print unknownListP[index].name
    #unknown positive ion mode
    k = len(dataListP)
    matchedNum = np.zeros(k)
    for mz in unknownListP[index].peak.sort_index(by="Relative Intensity", ascending=False)["m/z"][0:40]:
        count = 0
        for ref in dataListP:    
            if mz in array(ref.peak.sort_index(by="Relative Intensity", ascending=False)["m/z"][0:40]):
                matchedNum[count] += 1
            count += 1    
            
    print matchedNum            
    print dataListP[matchedNum.argmax()].name
    print "score: " ,  max(matchedNum)/numPeaks * 100
    
    #unknown negative ion mode
    k = len(dataListN)
    matchedNum = np.zeros(k)
    for mz in unknownListN[index].peak.sort_index(by="Relative Intensity", ascending=False)["m/z"][0:40]:
        count = 0
        for ref in dataListN:
        
            if mz in array(ref.peak.sort_index(by="Relative Intensity", ascending=False)["m/z"][0:40]):
                matchedNum[count] += 1
                
            count += 1    
    
    print matchedNum
    print dataListN[matchedNum.argmax()].name
    print "score: " ,  max(matchedNum)/numPeaks * 100
    print "\n"
    
    
for index in range(4):
    UnknownMatching(index)    



"""""
Blood sample matching
"""""
def BloodMatching(index):
    print bloodP[index].name
    
    #blood sample positive ion mode
    k = len(dataListP)
    matchedNum = np.zeros(k)
    for mz in bloodP[index].peak.sort_index(by="Relative Intensity", ascending=False)["m/z"][0:40]:
        count = 0
        for ref in dataListP:
        
            if mz in array(ref.peak.sort_index(by="Relative Intensity", ascending=False)["m/z"][0:40]):
                matchedNum[count] += 1
            count += 1    
            
    print matchedNum         
    print dataListP[matchedNum.argmax()].name
    print "score: " ,  max(matchedNum)/numPeaks * 100
    
    #blood sample negative ion mode
    k = len(dataListN)
    matchedNum = np.zeros(k)
    for mz in bloodN[index].peak.sort_index(by="Relative Intensity", ascending=False)["m/z"][0:40]:
        count = 0
        for ref in dataListN:
        
            if mz in array(ref.peak.sort_index(by="Relative Intensity", ascending=False)["m/z"][0:40]):
                matchedNum[count] += 1
            count += 1    
    
    print matchedNum
    print dataListN[matchedNum.argmax()].name
    print "score: " ,  max(matchedNum)/numPeaks * 100
    print "\n"

for index in range(4):
        BloodMatching(index)

"""
Self matching for debugging
"""

def SelfMatching(index):
    k = len(dataListP)
    matchedNum = np.zeros(k)
    for mz in dataListP[index].peak.sort_index(by="Relative Intensity", ascending=False)["m/z"][0:40]:
        count = 0
        for ref in dataListP:
            if mz in array(ref.peak.sort_index(by="Relative Intensity", ascending=False)["m/z"][0:40]):
                matchedNum[count] += 1
            count += 1    

    print matchedNum
    print dataListP[matchedNum.argmax()].name
    print "score: " ,  max(matchedNum)/numPeaks * 100
    
    k = len(dataListN)
    matchedNum = np.zeros(k)
    for mz in dataListN[index].peak.sort_index(by="Relative Intensity", ascending=False)["m/z"][0:15]:
        count = 0
        for ref in dataListN:
            if mz in array(ref.peak.sort_index(by="Relative Intensity", ascending=False)["m/z"][0:15]):
                matchedNum[count] += 1
            count += 1    

    print matchedNum
    print dataListN[matchedNum.argmax()].name
    print "score: " ,  max(matchedNum)/15 * 100
    print "\n"

for index in range(9):
    SelfMatching(index)





