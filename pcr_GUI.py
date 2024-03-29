# coding: utf-8
# PCR_recipe.pyと同じディレクトリにこのファイルを置く
from  main import pcr_recipe
import PySimpleGUI as sg
import pandas as pd
import datetime
import os
import shutil
from bs4 import BeautifulSoup
import webbrowser

def bind_html(html_path_list: list) -> str:
    '''複数のhtmlを順番に結合する

    Args:
        html_path_list ([str]): htmlファイルのpathリスト

    Attributes:
        pure_bound_html (str): 純粋にhtmlを文字列として 結合したもの
        bound_html (str): pure_bound_htmlから結合部分を 取り除いたもの

    Returns:
        str
    '''
    soup_list = []
    for html_path in html_path_list:
        with open(html_path, encoding="utf-8") as f:
            soup_list.append(BeautifulSoup(f.read()))
    pure_bound_html = ''.join([soup.prettify() for soup in soup_list])
    bound_html = pure_bound_html.replace('<table border="1" class="dataframe">', '<br>\n<table>')
    bound_html = bound_html.replace('<tr style="text-align: right;">', '<tr>')
    return bound_html

# レイアウト設計
layout_input =[
    [sg.Text("プライマー foward (5' → 3')",size=(25, 1)), sg.InputText(key = 'primer_fw')],
    [sg.Text("プライマー reverse (5' → 3')",size=(25, 1)), sg.InputText(key = 'primer_rv')],
    [sg.Text("プライマー濃度 (μM)",size=(25, 1)), sg.InputText("10", key = 'primer_conc_μM')],
    [sg.Text("増幅する領域",size=(25, 1)), sg.InputText(key = 'amplify_region')],
    [sg.Text("テンプレート濃度 (ng/μL)",size=(25, 1)), sg.InputText("1", key = 'template_conc_ng_μL')],
    [sg.Text("使うメーカー",size=(25, 1)), sg.Combo(["KOD -Plus-", "KOD One", "PrimeSTAR", "Ex Taq"])],
    [sg.Text("サンプルの本数",size=(25, 1)), sg.InputText("1", key = 'sample_size')],
    [sg.Text("反応総量 (μL)",size=(25, 1)), sg.InputText("25", key = 'total_vol_μL_per_sample')],
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

column_1 = [
    [sg.Text("プライマー forward: " + str(pcr.primer_fw_seq), size=(50, 2))],
    [sg.Text("プライマー reverse: " + str(pcr.primer_rv_seq), size=(50, 2))],
    [sg.Text("Tm値 (Wallace法): " + str(round(pcr.tm_value_Wallace, 1)) + "°C",size=(25, 1))],
    [sg.Text("Tm値 (GC法): " + str(round(pcr.tm_value_GC, 1)) + "°C",size=(25, 1))],
    [sg.Text("Tm値 (最近接塩基法, 推奨値): " + str(round(pcr.tm_value_NN, 1)) + "°C",size=(25, 1))],
    [sg.Text("試薬: " + str(pcr.reagent_name))]
    ]

column_2 = [
    [sg.Text("Termal cycler設定")],
    [sg.Table(
    values = pcr.thermal_list,
    headings = pcr.thermal_col_name)]
    ]

layout_tables = [
    [sg.Column(column_1, justification= "center"), sg.Column(column_2, justification= "center")],
    [sg.Table(
        values = pcr.conc_table_list,
        headings = pcr.conc_col_name)],
    [sg.Button("HTMLに出力"), sg.Button("保存せずに終了")]
    ]

window = sg.Window('PCR レシピ 結果', layout_tables, resizable = True)

html_script = \
    '<html>\n'\
    + '<head>\n' \
        + '<link rel="stylesheet" type="text/css" href="table_style.css">\n'\
    + '</head>\n'\
    + '<body>\n'\
    + '<br>\n'\
    + "プライマー forward: " + str(pcr.primer_fw_seq) \
    + "\n<br>"\
    + "\nプライマー reverse: " + str(pcr.primer_rv_seq) \
    + "\n<br>\n"\
    + "使用した酵素: " + str(pcr.reagent_name)

html_script_end  = \
    '\n</body>\n'\
    + '</html>'
# Event Loop to process "events" and get the "values" of the inputs
while True:
    event, values = window.read()
    if event in (sg.WIN_CLOSED, '保存せずに終了'):
        break
    if event in ("HTMLに出力"):
        os.mkdir("temp") # 一時フォルダの作成
        pcr.tm_table.to_html("temp/" + str(datetime.date.today()) + "_tm_table.html", index = False)
        pcr.thermal_table.to_html("temp/" + str(datetime.date.today()) + "_thermal_cycle.html", index = False)
        pcr.conc_table.to_html("temp/" + str(datetime.date.today()) + "_conc_table.html", index = False)
        html_path_list = ["temp/" + str(datetime.date.today()) + "_tm_table.html", "temp/" + str(datetime.date.today()) + "_thermal_cycle.html", "temp/" + str(datetime.date.today()) + "_conc_table.html"]
        bound_html = bind_html(html_path_list) # html結合
        with open(str(datetime.date.today())+".html", mode='w', encoding="utf-8") as f:
            f.write(html_script)
            f.write(bound_html) # html出力
            f.write(html_script_end)
        shutil.rmtree("temp/") # 一時フォルダの削除
        print(sg.popup('HTMLファイルを出力しました！\n' +"保存先: " + os.path.abspath(str(datetime.date.today()) + ".html")))
        html_path = os.path.abspath(str(datetime.date.today()) + ".html")
        break
window.close()