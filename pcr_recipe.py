# coding: utf-8
from Bio import SeqIO
import dataclasses
from dataclasses import dataclass
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import pandas as pd
import numpy as np

@dataclass #required python >= 3.7
class pcr_recipe: #インスタンス
    primer_fw: str = "" # primerの塩基配列
    primer_rv: str = ""
    amplify_region: str = "" # 増幅させる予定の塩基配列
    reagent_name: str = "" # "KOD" or "PrimeSTAR"
    template_conc_ng_μL: int = 1 # template DNAの濃度 1ng/μLが初期
    primer_conc_μM: int = 10 # primer DNAの濃度 10 μMが初期
    total_vol_μL_per_sample: int = 25 # 1サンプルあたりの反応液総量 25 μLが初期
    sample_size: int = 5 # サンプルサイズ
    
    def calc_tm_value(self):
        if self.primer_fw == 0 or self.primer_rv == 0:
            print("error！primerが入力されていません")
        else:
            self.primer_fw_seq = Seq(self.primer_fw.upper(), IUPAC.ambiguous_dna)
            self.primer_rv_seq = Seq(self.primer_rv.upper(), IUPAC.ambiguous_dna)
            self.tm_value_NN = np.mean([mt.Tm_NN(self.primer_fw_seq, Na = 50, nn_table = mt.DNA_NN1), mt.Tm_NN(self.primer_rv_seq, Na = 50, nn_table = mt.DNA_NN1)]) # 最近接塩基法で計算
            self.tm_value_GC = np.mean([mt.Tm_GC(self.primer_fw_seq, Na = 50, valueset = 1), mt.Tm_GC(self.primer_rv_seq, Na = 50, valueset = 1)]) # GC法で計算
    
    def count_amplify_region(self):
        self.amplify_region_seq = Seq(self.amplify_region.upper(), IUPAC.ambiguous_dna)
        self.amplify_regipon_kbp = len(self.amplify_region_seq)/1000
    
    def create_thermal_cycle_plan(self):
        self.calc_tm_value()
        self.count_amplify_region()
        if self.reagent_name == "KOD": # KODの時の温度設定
            self.thermal_list = [
                ["94°C", "2 min"],
                ["94°C", "15 sec"],
                [str(round(self.tm_value_NN - 5)) + "°C", "30 sec"],
                ["68°C",  str(round(self.amplify_regipon_kbp, 2)) + " min"]
                ]
                
        elif self.reagent_name == "PrimeSTAR": # PrimeSTARの時の温度設定
            if self.tm_value_GC - 5 >= 55:
                self.thermal_list = [
                ["98°C", "10 sec"],
                ["55°C", "5 sec"],
                ["72°C", str(round(self.amplify_regipon_kbp, 2)*5) + "sec"],
                ]
            else:
                self.thermal_list = [
                ["98°C", "10 sec"],
                ["55°C", "15 sec"],
                ["72°C", str(round(self.amplify_regipon_kbp, 2)*5) + "sec"],
                ]
                self.message_tm =print(
                    "PCR thermal setting\n" + "98°C: 10 sec\n" + "55℃: 15 sec\n"
                    + "72°C: " + str(round(self.amplify_regipon_kbp, 2)*5) + " sec"
                    )
        self.thermal_col_name = ["PCR temperature","time"]
        self.thermal_table= pd.DataFrame(self.thermal_list, columns = self.thermal_col_name)
    
    def create_pcr_recipe(self): # 雑に作りすぎた。要改善
        self.create_thermal_cycle_plan()
        self.conc_col_name = [r"Reagent", r"Usage/sample (μL)", r"必要量 (μL)", "Final concentration"]
        self.sample_amount = self.total_vol_μL_per_sample*self.sample_size*1.2
        if self.reagent_name == "KOD":
            self.conc_list = [
                ["10×PCR Buffer", self.total_vol_μL_per_sample/10, self.sample_amount/10, "1x"],
                ["25 mM MgSO4 solution", round(self.total_vol_μL_per_sample/16.7, 2), round(self.sample_amount/16.7, 2), "1.5 mM"],
                [str(self.primer_conc_μM) + " μM Primer_fw", self.total_vol_μL_per_sample/(self.primer_conc_μM/0.3), round(self.sample_amount/(self.primer_conc_μM/0.3), 3), "0.3 μM"],
                [str(self.primer_conc_μM) + " μM Primer_rv", self.total_vol_μL_per_sample/(self.primer_conc_μM/0.3), round(self.sample_amount/(self.primer_conc_μM/0.3), 3), "0.3 μM"],
                ["2 mM dNTPs", self.total_vol_μL_per_sample/10, self.sample_amount/10, "0.2 mM"],
                [str(self.template_conc_ng_μL) + " ng/μL Template DNA", self.total_vol_μL_per_sample/(self.template_conc_ng_μL/0.04), self.sample_amount/(self.template_conc_ng_μL/0.04), "0.04 ng/μL"],
                ["KOD plus (1U/μL)", self.total_vol_μL_per_sample/50, self.sample_amount/50, "1 U"]
            ]
        elif self.reagent_name == "PrimeSTAR":
            self.conc_list = [
                ["PrimeSTAR Max Premix (2×)", self.total_vol_μL_per_sample/2, self.sample_amount/2, "1x"],
                [str(self.primer_conc_μM) + " μM Primer_fw", self.total_vol_μL_per_sample/(self.primer_conc_μM/0.3), round(self.sample_amount/(self.primer_conc_μM/0.3), 3), "0.3 μM"],
                [str(self.primer_conc_μM) + " μM Primer_rv", self.total_vol_μL_per_sample/(self.primer_conc_μM/0.3), round(self.sample_amount/(self.primer_conc_μM/0.3), 3), "0.3 μM"],
                [str(self.template_conc_ng_μL) + " ng/μL Template DNA", self.total_vol_μL_per_sample/(self.template_conc_ng_μL/0.04), self.sample_amount/(self.template_conc_ng_μL/0.04), "0.04 ng/μL"],
            ]
        self.conc_table= pd.DataFrame(self.conc_list, columns= self.conc_col_name)
        self.dw_list = ["DW", self.total_vol_μL_per_sample - self.conc_table["Usage/sample (μL)"].sum(), self.sample_amount - self.conc_table["必要量 (μL)"].sum(), "-"]
        self.conc_table_total = ["Total", self.total_vol_μL_per_sample,  self.sample_amount, "-"]
        self.conc_table.loc[7] = self.dw_list
        self.conc_table.loc[8] = self.conc_table_total
        self.conc_table['sample size'] = self.sample_size
        self.conc_table_list = self.conc_table.values.tolist()