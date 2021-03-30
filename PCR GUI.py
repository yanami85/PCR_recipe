# coding: utf-8
# PCR_recipe.pyと同じディレクトリにこのファイルを置く
from  pcr_recipe import pcr_recipe
import PySimpleGUI as sg
import matplotlib.pyplot as plt
import pandas as pd
import datetime

# レイアウト設計
layout_input =[
    [sg.Text("プライマー foward (5' → 3')",size=(20,1)), sg.InputText('ttaatacgactcactatagg', key = 'primer_fw')],
    [sg.Text("プライマー reverse (5' → 3')",size=(20,1)), sg.InputText('ttaatacgactcactatagg', key = 'primer_rv')],
    [sg.Text("プライマー濃度 (μM)",size=(20,1)), sg.InputText("10", key = 'primer_conc_μM')],
    [sg.Text("増幅する領域",size=(20,1)), sg.InputText("A"*10000, key = 'amplify_region')],
    [sg.Text("テンプレート濃度 (ng/μL)",size=(20,1)), sg.InputText("1", key = 'template_conc_ng_μL')],
    [sg.Text("使うメーカー",size=(20,1)), sg.Combo(["KOD -Plus-", "KOD One", "PrimeSTAR"])],
    [sg.Text("サンプルの本数",size=(20,1)), sg.InputText("1", key = 'sample_size')],
    [sg.Text("反応総量 (μL)",size=(20,1)), sg.InputText("25", key = 'total_vol_μL_per_sample')],
    [sg.Submit(button_text = "Cancel"), sg.Submit(button_text = "次へ")]
    ]
# 窓を作る
window = sg.Window('PCR レシピ', layout_input, resizable = True)

# Event Loop to process "events" and get the "values" of the inputs
while True:
    event, values = window.read()
    if event in (sg.WIN_CLOSED, 'Cancel'):
        exit()
    if event == "次へ":
        break
window.close()

# import pcr_recipe内で定義した変数を入力
pcr = pcr_recipe()
pcr.primer_fw = values["primer_fw"]
pcr.primer_rv = values["primer_rv"]
pcr.amplify_region = values["amplify_region"]
pcr.reagent_name = values[0]
pcr.template_conc_ng_μL = int(values["template_conc_ng_μL"])
pcr.primer_conc_μM = int(values["primer_conc_μM"])
pcr.sample_size = int(values["sample_size"])
pcr.total_vol_μL_per_sample = int(values["total_vol_μL_per_sample"])

# import pcr_recipe内で定義した関数を実行
pcr.create_pcr_recipe()

layout_tables = [
    [sg.Text("Tm値 (for KOD): " + str(round(pcr.tm_value_NN, 1)) + "°C",size=(30,1))],
    [sg.Text("Tm値 (for PrimeSTAR): " + str(round(pcr.tm_value_GC, 1)) + "°C",size=(30,1))],
    [sg.Table(
        key='-TABLE-',
        values = pcr.thermal_list,
        headings = pcr.thermal_col_name
        )],
    [sg.Text(str(pcr.reagent_name),size=(20,1))],
    [sg.Table(
        key='-TABLE-',
        values = pcr.conc_table_list,
        headings = pcr.conc_col_name
        )],
    [sg.Button("Exit"), sg.Button("HTMLに出力")],
    ]

window = sg.Window('pcr_recipe', layout_tables, resizable = True)


# Event Loop to process "events" and get the "values" of the inputs
while True:
    event, values = window.read()
    if event in (sg.WIN_CLOSED, 'Exit'):
        break
    if event in ("HTMLに出力"):
        pcr.thermal_table.to_html(str(datetime.date.today()) + "_thermal_cycle.html")
        pcr.conc_table.to_html(str(datetime.date.today()) + "_conc_table.html")
        popup = sg.popup_ok('HTMLファイルを出力しました！')
        print(popup)
window.close()
