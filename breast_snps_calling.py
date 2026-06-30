
import subprocess

inp1=input('Whether rs2366152 SNP exists? (Y/N)')
while(1):
    if(inp1=='Y'):
        subprocess.run(['python', 'breast/rs2366152_predict.py'])
    elif(inp1=='N'):
        print('rs2366152 SNP does not exists.')
    else:
        print('Only specify the mentioned Symbol.')
    
    if ((inp1=='Y')|(inp1=='N')):
        break
    else:
        inp1=input('Whether rs2366152 SNP exists? (Y/N)')
    

print('#####################################\n')

inp2=input('Whether rs7091441 SNP exists? (Y/N)')
while(1):
    if(inp2=='Y'):
        subprocess.run(['python', 'breast/rs7091441_predict.py'])
    elif(inp2=='N'):
        print('rs7091441 SNP does not exists.')
    else:
        print('Only specify the mentioned Symbol.')
    
    if ((inp2=='Y')|(inp2=='N')):
        break
    else:
        inp2=input('Whether rs7091441 SNP exists? (Y/N)')
    


