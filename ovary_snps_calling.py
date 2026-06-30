
import subprocess

inp1=input('Whether rs932501 SNP exists? (Y/N)')
while(1):
    if(inp1=='Y'):
        subprocess.run(['python', 'ovary/rs932501_predict.py'])
    elif(inp1=='N'):
        print('rs932501 SNP does not exists.')
    else:
        print('Only specify the mentioned Symbol.')
    
    if ((inp1=='Y')|(inp1=='N')):
        break
    else:
        inp1=input('Whether rs932501 SNP exists? (Y/N)')
    

print('#####################################\n')

inp2=input('Whether rs2072588 SNP exists? (Y/N)')
while(1):
    if(inp2=='Y'):
        subprocess.run(['python', 'ovary/rs2072588_predict.py'])
    elif(inp2=='N'):
        print('rs2072588 SNP does not exists.')
    else:
        print('Only specify the mentioned Symbol.')
    
    if ((inp2=='Y')|(inp2=='N')):
        break
    else:
        inp2=input('Whether rs2072588 SNP exists? (Y/N)')
    

print('#####################################\n')

inp3=input('Whether rs7318592 SNP exists? (Y/N)')
while(1):
    if(inp3=='Y'):
        subprocess.run(['python', 'ovary/rs7318592_predict.py'])
    elif(inp3=='N'):
        print('rs7318592 SNP does not exists.')
    else:
        print('Only specify the mentioned Symbol.')
    
    if ((inp3=='Y')|(inp3=='N')):
        break
    else:
        inp3=input('Whether rs7318592 SNP exists? (Y/N)')

print('#####################################\n')

inp4=input('Whether rs9506960 SNP exists? (Y/N)')
while(1):
    if(inp4=='Y'):
        subprocess.run(['python', 'ovary/rs9506960_predict.py'])
    elif(inp4=='N'):
        print('rs9506960 SNP does not exists.')
    else:
        print('Only specify the mentioned Symbol.')
    
    if ((inp4=='Y')|(inp4=='N')):
        break
    else:
        inp4=input('Whether rs9506960 SNP exists? (Y/N)')


print('#####################################\n')

inp5=input('Whether rs9510420 SNP exists? (Y/N)')
while(1):
    if(inp5=='Y'):
        subprocess.run(['python', 'ovary/rs9510420_predict.py'])
    elif(inp5=='N'):
        print('rs9510420 SNP does not exists.')
    else:
        print('Only specify the mentioned Symbol.')
    
    if ((inp5=='Y')|(inp5=='N')):
        break
    else:
        inp5=input('Whether rs9510420 SNP exists? (Y/N)')

print('#####################################\n')

inp6=input('Whether rs12583808 SNP exists? (Y/N)')
while(1):
    if(inp6=='Y'):
        subprocess.run(['python', 'ovary/rs12583808_predict.py'])
    elif(inp6=='N'):
        print('rs12583808 SNP does not exists.')
    else:
        print('Only specify the mentioned Symbol.')
    
    if ((inp6=='Y')|(inp6=='N')):
        break
    else:
        inp6=input('Whether rs12583808 SNP exists? (Y/N)')

print('#####################################\n')

inp7=input('Whether rs60135126 SNP exists? (Y/N)')
while(1):
    if(inp7=='Y'):
        subprocess.run(['python', 'ovary/rs60135126_predict.py'])
    elif(inp7=='N'):
        print('rs60135126 SNP does not exists.')
    else:
        print('Only specify the mentioned Symbol.')
    
    if ((inp7=='Y')|(inp7=='N')):
        break
    else:
        inp7=input('Whether rs60135126 SNP exists? (Y/N)')
