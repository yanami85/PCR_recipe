# PCR_recipe

## What's PCR_recipe?
PCRで必要な各種試薬の量やthermal cyclerの温度設定を自動的に出力するソフト。
**(python > 3.7)**

## 必要なもの
+ dataclasses
+ Biopython
+ pandas
+ numpy
+ pySimpleGUI

pip or anacondaで適宜入れて下さい。
## 使用方法
1. PCR GUI.pyとpcr_recipe.pyを同じディレクトリに格納する。
2. PCR GUI.pyを走らせる
3. 以下の画面が出てくるので、primer、プライマー濃度、増幅領域、テンプレート濃度 etc. を入力
(使うメーカーは、2020-03-29現在で**KOD Plus, KOD One, PrimeSTAR**に対応している。) 
![image](https://user-images.githubusercontent.com/41857834/113079382-e2131c00-920f-11eb-82f9-49685069df47.png)
4. 「次へ」をクリックすれば、Thermal cyclerの設定と必要試薬量が一覧で見られる。
![image](https://user-images.githubusercontent.com/41857834/113079403-eb9c8400-920f-11eb-9aa4-2a7bbfac21dc.png)
5. 「HTMLに出力」をクリックすると1. のディレクトリにHTMLファイルが保存される
![image](https://user-images.githubusercontent.com/41857834/113079291-bbed7c00-920f-11eb-86d8-ef2df6d347d5.png)
