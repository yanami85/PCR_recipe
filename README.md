# PCR_recipe

## What's PCR_recipe?
PCRで必要な各種試薬の量やthermal cyclerの温度設定を自動的に出力するソフト。
**python 3.7 or 3.8を推奨します。**

## 必要なもの
+ dataclasses 0.6
+ Biopython 1.78
+ pandas 1.2.3
+ numpy 1.19.2
+ pySimpleGUI 4.38.0

pip or anacondaで適宜入れて下さい。
(**python 3.9だとBiopython 1.78をインストールできない可能性があります。**)

## 使用方法
1. PCR GUI.pyとpcr_recipe.py, df_style.cssを同じディレクトリに格納する。
2. PCR GUI.pyを走らせる
3. 以下の画面が出てくるので、primer、プライマー濃度、増幅領域、テンプレート濃度 etc. を入力
(使うメーカーは、2020-03-29現在で**KOD Plus, KOD One, PrimeSTAR**に対応している。) 

![image](https://user-images.githubusercontent.com/41857834/113080175-6ade8780-9211-11eb-9497-d704d9e512cd.png)

4. 「次へ」をクリックすれば、Thermal cyclerの設定と必要試薬量が一覧で見られる。

![image](https://user-images.githubusercontent.com/41857834/113080183-729e2c00-9211-11eb-9556-806c78aa8e35.png)

5. 「HTMLに出力」をクリックすると1. のディレクトリにHTMLファイルが保存される

![image](https://user-images.githubusercontent.com/41857834/113080229-85186580-9211-11eb-9754-1793dae8d3b2.png)
