
str1=['D:\lib\site-packages\pymoo\\algorithms\\nsga2.py',
'D:\lib\site-packages\pymoo\\algorithms\genetic_algorithm.py',
'D:\lib\site-packages\pymoo\model\mating.py',
'D:\lib\site-packages\pymoo\model\infill.py']

str2=['nsga2.py','genetic_algorithm.py','mating.py','infill.py']
for i in range(len(str1)):
    fp = open(str1[i], "w")
    fp.truncate()
    with open(str2[i], "r") as fileobj:
        text = fileobj.read()
    with open(str1[i], "w") as fw:
        fw.write(text)

