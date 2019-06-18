import PySimpleGUI as sg

def GUI():

    #I know that this is written terribly, its my first try to with a GUI

    startupLayout = [[sg.Text('Enter the number gating conditions you want:')],
                     [sg.Text('(enter 0 if you just want the raw data)')],
                     [sg.Input(do_not_clear = True,key='input'),sg.Submit(),sg.Cancel()]]

    startup = sg.Window('Gate number selection').Layout(startupLayout)
    eventS,valuesS = startup.Read()
    gatNum = int(valuesS['input'])
    startup.Close()

    if gatNum == 0:

        layout = [[sg.Text('Outputting raw data to csv...')],
                  [sg.Text('Select the file to open:'),sg.FileBrowse(key='fileIn')],
                  [sg.Text('Select output folder:'),sg.FolderBrowse(key='foldOut')],
                  [sg.Text('Write output file name (end with .csv):'),sg.Input(key='nameOut')],
                  [sg.Submit(), sg.Cancel()]]              
        

    elif isinstance(gatNum, int):

        layout = [[sg.Text("Outputting %s csv(s)."% gatNum)],      
                  [sg.Text('Select the file to open:'),sg.FileBrowse(key='fileIn')],
                  [sg.Text('Select output folder:'),sg.FolderBrowse(key='foldOut')],
                  [sg.Text('Write output file name (end with .csv):'),sg.Input(key='nameOut')],
                  [sg.Text('')], 
                  [sg.Text('Enter Gate(s) as [difmin,difmax,summin,summax]')],
                  *[[sg.Input(do_not_clear = True,key='bound%s'%i)] for i in range(gatNum)],
                                    [sg.Submit(), sg.Cancel()]]      


    window=sg.Window('Choose Raw File').Layout(layout)
    event,values = window.Read()
    window.Close()
    
    fileIn = values['fileIn']
   
    fileOut= values['foldOut']+'/'+values['nameOut']

    if gatNum != 0:
        channelBounds = [ v for k,v in values.items() if k.startswith('bound')]
    else:
        channelBounds = None

    return fileIn,fileOut,channelBounds 

