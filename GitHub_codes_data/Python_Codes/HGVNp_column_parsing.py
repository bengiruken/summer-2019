#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 18:07:00 2019

@author: bengi
"""

'''HGVNp_short column parsing '''
#Sadece point mutation olan line'lar ve mutasyondan sonra oluşan residue'nun ne olduğuna dair bilgi varsa
#line'ı yazdırıyoruz
'''MUTATION CASES
*First two cases are the same 
1.)Frame_Shift_Del
HGVNp_short_parse("p.G80Vfs*3")
('80', 'G', 'Vfs*3')[2].split('fs*')[0]
2.)Frame_Shift_Ins
HGVNp_short_parse("p.F439Ifs*2")
3.)In_Frame_Del: _: indicates a residue interval
    *exlude this mutation type
    p.V23249_T23250del
4.)In_Frame_Ins
    *effects more than one residue position, so exclude???
    p.F32delinsCI
    p.Q329_Q330dup
5.)Intron
    HGVS is just ".", exclude
6.)Missense mutation, include
7.)Nonsense_Mutation, include
    p.E956* : here * is the termination codon, write a case for this
    p.Y216_K218delins* : should not include this, exclude containing "_" ones
8.)Nonstop_Mutation: aşağıdaki case'e göre al
    p.*770Yext*28
 HGVNp_short_parse("p.*770Yext*28")
('770', '*', 'Yext*28')[2].split('ext*')[0]
9.)RNA, exclude
10)Splice_Site, exclude
    p.X437_splice
11)Translation_Start_Site, exclude
    p.M1?
'''

'''parsing hgvnp column'''
def HGVNp_short_parse(hgvn):
    import re
    test_string=hgvn
    a=[int(s) for s in re.findall(r'-?\d+\.?\d*', test_string)] 
    if len(a)<1:
        return None,None,None
    else:
        a=str(a[0])
        u=test_string.find(a[0])
        v=2+len(a)
        aminoacid=test_string[2:u]
        mutant=test_string[v+1:]
    return a,aminoacid,mutant  
HGVNp_short_parse("p.*770Yext*28")
('770', '*', 'Yext*28')[2].split('ext*')[0]
HGVNp_short_parse("p.F439Ifs*2")
('439', 'F', 'Ifs*2')[2].split('fs*')[0]
