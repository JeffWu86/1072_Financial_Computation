HW2 
執行方法： ./HW2_1 input_1.txt
編譯: g++ HW2_1.cpp -o HW2_1

HW2_1: All
HW2_2: Change (n m)p^n-j*(1-p)^j
HW2_3: 改成vector版本 可以讓ｎ跑到3000 且搭配binomial可以算出結果 if利用第一版本算Cn取j會爆掉
n=30000逼近BS答案到小數點第三位
HW2_4: 刪除basic bonus1 避免跑太久

Input.txt:
S K r q sigma T #_of_sim #_of_rep n