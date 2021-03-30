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

![image](https://user-images.githubusercontent.com/41857834/112975935-3629ec00-918f-11eb-8fc1-2d04f53a6c8c.png)

4. 「次へ」をクリックすれば、Thermal cyclerの設定と必要試薬量が一覧で見られる。
![image](https://user-images.githubusercontent.com/41857834/112975970-3f1abd80-918f-11eb-9953-96c552083616.png)

5. 「HTMLに出力」をクリックすると1. のディレクトリにHTMLファイルが保存される。
![image](https://user-images.githubusercontent.com/41857834/112976008-4a6de900-918f-11eb-8cc8-a16446bc6251.png)
