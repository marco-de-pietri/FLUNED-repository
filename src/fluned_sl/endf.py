"""
class to operate with ENDF files - authors j. Alguacil UNED, P. Sauvan UNED
"""
from math import log10

def z_sym(Z):
    """
    function that returns the atomic numbe of a atomic symbol
    """
    dict_z=({ 1:'H' , 2:'He', 3:'Li', 4:'Be', 5:'B' , 6:'C' , 7:'N' , 8:'O' ,
           	  9:'F' ,10:'Ne',11:'Na',12:'Mg',13:'Al',14:'Si',15:'P' ,16:'S' ,
             17:'Cl',18:'Ar',19:'K' ,20:'Ca',21:'Sc',22:'Ti',23:'V' ,24:'Cr',
             25:'Mn',26:'Fe',27:'Co',28:'Ni',29:'Cu',30:'Zn',31:'Ga',32:'Ge',
             33:'As',34:'Se',35:'Br',36:'Kr',37:'Rb',38:'Sr',39:'Y' ,40:'Zr',
             41:'Nb',42:'Mo',43:'Tc',44:'Ru',45:'Rh',46:'Pd',47:'Ag',48:'Cd',
             49:'In',50:'Sn',51:'Sb',52:'Te',53:'I' ,54:'Xe',55:'Cs',56:'Ba',
             57:'La',58:'Ce',59:'Pr',60:'Nd',61:'Pm',62:'Sm',63:'Eu',64:'Gd',
             65:'Tb',66:'Dy',67:'Ho',68:'Er',69:'Tm',70:'Yb',71:'Lu',72:'Hf',
             73:'Ta',74:'W' ,75:'Re',76:'Os',77:'Ir',78:'Pt',79:'Au',80:'Hg',
             81:'Tl',82:'Pb',83:'Bi',84:'Po',85:'At',86:'Rn',87:'Fr',88:'Ra',
             89:'Ac',90:'Th',91:'Pa',92:'U' ,93:'Np',94:'Pu',95:'Am',96:'Cm',
             97:'Bk',98:'Cf',99:'Es',100:'Fm',101:'Md',102:'No'})

    return dict_z[Z]

def endf_format(x):
    """
    function to format correctly numbers for endf files
    """
    if  isinstance(x,int)  :
        return f'{x:11d}'
    if x == 0 :
        return " 0.000000+0"
    if x == "blnk" :
        return "           "
    if  isinstance(x,str)  :
        try :
            x=float(x)
        except TypeError:
            return x

        exp=log10(abs(x))
        expi=int(exp)
        if (exp < 0) :
            expi=expi-1
        flt=x*10**(-expi)
        if (flt >= 10.):
            flt=flt/10
            expi=expi+1
        if exp < -9 :
            return f'{flt:8.5f}{expi:3d}'
        elif -9 <= exp < 0 :
            return f'{flt:9.6f}{expi:2d}'
        elif  0 <= exp < 10 :
            return f'{flt:9.6f}+{expi:d}'
        elif 10 <= exp  :
            return f'{flt:8.5f}+{expi:2d}'

def get_value(x):
    """
    function description
    """

    if isinstance(x,int):
        return x
    if isinstance(x,float):
        return x
    if x.strip() == '':
        return 'blnk'

    try:
        x=int(x)
        return x
    except:
        try:
            x=float(x)
            return x
        except:
            if  x[-2] == '+' or x[-2]=='-':
                if(x[-1]==' ') :
                    exp=x[-2]+'0'
                else:
                    exp=x[-2:]
                xt=x[0:-2]+'E'+exp
            else:
                if(x[-2]==' ') :
                    exp=x[-3]+'0'+x[-1]
                else:
                    exp=x[-3:]
                xt=x[0:-3]+'E'+exp
            return float(xt)


def sum_tab1(xstab):
    """
    function description
    """

    np=len(xstab[0])
    sum_val=[0.]*np
    for xs in xstab:
        for i,val in enumerate(xs):
            sum_val[i]=sum_val[i]+val[1]

    xstot=[]
    for i,val in enumerate(xstab[0]):
        xstot.append([val[0],sum_val[i]])

    return xstot

def sort_MFMT_old(tab):
    """
    function description
    """

    sorted_list=[]
    for key in tab:
        if 'end' in key:
            continue
        if key[3] == '_' :
            k=3
        else:
            k=4
        MF=int(key[2:k])
        MT=int(key[k+3:])
        sorted_list.append(MF*1000+MT)
    sorted_list.sort()
    keys=[]
    for val in sorted_list:
        MF=val/1000
        MT=val%1000
        MFMT=f'MF{MF:d}_MT{MT:d}' %(MF,MT)
        keys.append([MFMT,MF,MT])
    return keys

def sort_MFMT(entry):
    """
    function description
    """
    sorted_list=[]
    invert={}
    for MF_MT in entry.keys():
        if 'end' not in MF_MT:
            ind=entry[MF_MT][0]*1000+entry[MF_MT][1]
            invert[ind]=MF_MT
        tab=invert.keys()

    tab.sort()
    for key in tab:
        MFMT=invert[key]
        sep=MFMT.index('_')
        MF=int(MFMT[2:sep])
        MT=int(MFMT[sep+3:])
        sorted_list.append([MFMT,MF,MT])
    return sorted

class endf_section:
    """
    class description
    """
    def __init__(self,ZA,AWR,MAT,MF,MT,IS):
        """
        function description
        """
        self.data=[]
        self.ZA=ZA
        self.AWR=AWR
        self.MAT=MAT
        self.MF=MF
        self.MT=MT
        self.IS=IS
        self.pointer=0

    def remove(self,pointer,nlines=1):
        """
        function description
        """
        del self.data[pointer:pointer+nlines]

    def insert(self,pointer,tab):
        """
        function description
        """
        self.data[pointer:pointer]=tab

    def insert_MT(self,MT,tab):
        """
        function description
        """
        self.IS=0
        self.MT=MT
        self.data=tab[:]
        nline=len(tab)
        try:
            self.fill_MF(MT,nline)
        except:
            return

    def fill_MF(self,MT,nline):
        """
        function description
        """
        if  len(self.entry.keys()) != 0 :
            MT_min=0
            iMT=-1
            for key in self.entry.keys():
                i1=key.index('_')+3
                if (key[i1:]=='d'):
                    key_MT=1000
                else:
                    key_MT=int(key[i1:])

                if  MT_min < key_MT < MT  :
                    iMT=self.entry[key][0]
                    MT_min=key_MT

            iMT=iMT+1
            MF_MT= f'MF{self.MF:d}_MT{MT:d}'
            self.filedata[iMT:iMT]=[self.data]

            for key in self.entry.keys():
                if (self.entry[key][0] >= iMT ):
                    self.entry[key][0] = self.entry[key][0]+1
            self.entry[MF_MT]=[iMT,nline]
        else:
            self.filedata.append(self.data[:])
            MF_MT= f'MF{self.MF:d}_MT{MT:d}'
            self.entry[MF_MT]=[0,nline]
        self.pos=self.pos+1


    def write_line(self,tab6,str_out=True):
        """
        function description
        """
        l6=''
        line=''
        for i,val in enumerate(tab6):
            l6=l6+endf_format(val)
            if ((i+1)%6 == 0):
                line=line+'{l6:s}{self.MAT:4d}{self.MF:2d}{self.MT:3d}{self.IS:5d}\n'
                l6=''
                self.IS=self.IS+1

        if ((i+1)%6 != 0):
            nb=6-((i+1)%6)
            for j in range(nb) :
                l6=l6+endf_format('blnk')

            line=line+'{l6:s}{self.MAT:4d}{self.MF:2d}{self.MT:3d}{self.IS:5d}\n'
            self.IS=self.IS+1

        if str_out:
            return line
        else:
            ltab=line.splitlines()
            tab=[]
            for line in ltab:
                tab.append(line+'\n')
            return tab

    def get_tab6(self,line):
        """
        function description
        """
        tab6=[]
        i1=0
        for i in range(6):
            i2=i1+11
            x_str=line[i1:i2]
            tab6.append(get_value(x_str))
            i1=i2

        return tab6

    def get_MAT(self,line):
        """
        function description
        """
        MAT=int(line[66:70])
        return MAT

    def get_MF(self,line):
        """
        function description
        """
        MF=int(line[70:72])
        return MF

    def get_MT(self,line):
        """
        function description
        """
        MT=int(line[72:75])
        return MT

    def get_MT_list(self,MTsort=True):
        """
        function description
        """
        MT_tab=[]

        if MTsort:
            for key in self.entry.keys():
                MT_str=key[key.index('_'):]
                if 'end' not in MT_str:
                    MT_tab.append(int(MT_str[3:]))
            MT_tab.sort()

        else:
            invert={}
            for MF_MT in self.entry.keys():
                MT_str=MF_MT[MF_MT.index('_'):]
                if 'end' not in MT_str:
                    invert[self.entry[MF_MT][1]]=int(MT_str[3:])
            tab=invert.keys()
            tab.sort()
            for key in tab:
                MT_tab.append(invert[key])

        return MT_tab

    def jump(self,nval,cont=False):
        """
        function description
        """
        nline= (nval-1)/6 + 1
        self.pointer=self.pointer+nline
        if (not cont) : self.pointer=self.pointer+1

    def get_list(self,i0=0):
        """
        function description
        """
        line=self.data[i0]
        tab6=self.get_tab6(line)
        C1=get_value(tab6[0])
        C2=get_value(tab6[1])
        L1=get_value(tab6[2])
        L2=get_value(tab6[3])
        nr= get_value(tab6[4])
        N2= get_value(tab6[5])

        nrline= (nr-1)/6 + 1
        lst=[]
        i0=i0+1
        for i in range(nrline):
            line=self.data[i0+i]
            tab6=self.get_tab6(line)
            for v in tab6:
                lst.append(v)
        return [C1,C2,L1,L2,nr,N2],lst,i0+nrline


    def get_tab1(self,i0=0):
        """
        function description
        """

        line=self.data[i0]
        tab6=self.get_tab6(line)
        C1=get_value(tab6[0])
        C2=get_value(tab6[1])
        L1=get_value(tab6[2])
        L2=get_value(tab6[3])
        nr= get_value(tab6[4])
        np= get_value(tab6[5])
        nrline= (nr-1)/3 + 1
        npline= (np-1)/3 + 1
        nbt=[]
        tab=[]
        tab1=[]
        i0=i0+1
        for i in range(nrline):
            line=self.data[i0+i]
            tab6=self.get_tab6(line)
            for v in tab6:
                tab1.append(v)
        for i in range(nr):
            nbt.append([tab1[2*i],tab1[2*i+1]])

        tab1=[]
        i0=i0+nrline
        for i in range(npline):
            line=self.data[i0+i]
            tab6=self.get_tab6(line)
            for v in tab6:
                tab1.append(v)
        for i in range(np):
            tab.append([tab1[2*i],tab1[2*i+1]])
        return [C1,C2,L1,L2,nr,np],[nbt,tab],i0+npline

    def get_tab2(self,i0=0):
        """
        function description
        """
        line=self.data[i0]
        tab6=self.get_tab6(line)
        C1=get_value(tab6[0])
        C2=get_value(tab6[1])
        L1=get_value(tab6[2])
        L2=get_value(tab6[3])
        nr= get_value(tab6[4])
        nz= get_value(tab6[5])

        nrline= (nr-1)/3 + 1
        nbt=[]
        tab1=[]
        i0=i0+1
        for i in range(nrline):
            line=self.data[i0+i]
            tab6=self.get_tab6(line)
            for v in tab6:
                tab1.append(v)
        for i in range(nr):
            nbt.append([tab1[2*i],tab1[2*i+1]])
        return [C1,C2,L1,L2,nr,nz],nbt,i0+nrline

    def title(self,string):
        """
        function description
        """
        line='{0:<66s}{1:>4d}{2:>2d}{3:>3d}{4:>5d}\n'.format(string,1,0,0,0)
        self.IS=1
        self.data.append(line)
        return

    def string(self,string):
        """
        function description
        """
        line='{0:<66s}{1:>4d}{2:>2d}{3:>3d}{4:>5d}\n'.format(string,self.MAT,
                                                             self.MF,self.MT,
                                                             self.IS)
        self.IS=self.IS+1
        self.data.append(line)
        return

    def head(self,tab6):
        """
        function description
        """

        headline=self.write_line(tab6)
        self.data.append(headline)
        return

    def send(self):
        """
        function description
        """

        sendline = f' 0.000000+0 0.000000+0          0          0          0          0{self.MAT:4d}{self.MF:2d}  099999\n'
        self.data.append(sendline)
        self.IS=0
        return

    def list(self,tab6,tab):
        """
        function description
        """

        tab6[4]=len(tab)
        line=self.write_line(tab6)
        self.data.append(line)

        if type(tab[0]) is not list :
            ltab=self.write_line(tab,str_out=False)
        else:
            tabxny=[]
            for E_tab in enumerate(tab):
                for v in E_tab:
                    tabxny.append(v)
            ltab=self.write_line(tabxny,str_out=False)
        for line in ltab:
            self.data.append(line)

        return


    def tab1(self,tab6,nbt,tab):
        """
        function description
        """

        tab6[4]=len(nbt)
        tab6[5]=len(tab)

        l0=self.write_line(tab6)
        self.data.append(l0)

        tabxy=[]
        for i,val in enumerate(nbt):
            tabxy.append(val[0])
            tabxy.append(val[1])

        ltab=self.write_line(tabxy,str_out=False)
        for line in ltab:
            self.data.append(line)

        tabxy=[]
        for i,val in enumerate(tab):
            tabxy.append(val[0])
            tabxy.append(val[1])

        ltab=self.write_line(tabxy,str_out=False)
        for line in ltab:
            self.data.append(line)

        return

    def tab2(self,tab6,nbt):
        """
        function description
        """

        tab6[4]=len(nbt)

        l0=self.write_line(tab6)
        self.data.append(l0)

        tabxy=[]
        for i,val in enumerate(nbt):
            tabxy.append(val[0])
            tabxy.append(val[1])

        ltab=self.write_line(tabxy,str_out=False)
        for line in ltab:
            self.data.append(line)

        return

class MF1(endf_section):
    """
    class description
    """

    def __init__(self,ZA=None,AWR=None,MAT=None):
        """
        function description
        """
        endf_section.__init__(self,ZA,AWR,MAT,1,None,0)
        self.filedata=[]
        self.ZA=ZA
        self.AWR=AWR
        self.MAT=MAT
        self.MF=1
        self.MT=451
        self.IS=0
        self.entry={}
        self.pos=0

    def set_MF(self,tab,entry):
        """
        function description
        """
        self.filedata=tab[:]
        self.entry=entry

    def MF_end(self):
        """
        function description
        """
        fendline= f' 0.000000+0 0.000000+0          0          0          0          0{self.MAT:4d} 0  0    0\n'
        self.filedata.append(fendline)
        MF_end= f'MF{self.MF:d}_end'
        self.entry[MF_end]=[self.pos,1]

    def build_info(self,title,infofile,content,iso=0):
        """
        function description
        """

        self.title(title)
        LRP=0
        AWI=1.   # incident particle mass (in neutron mass)
        LISO=iso   #isomeric state

        NWD=len(infofile)
        NXC=len(content)

        tab6=[self.ZA,self.AWR,LRP,0,0,1]
        self.head(tab6)

        tab6=[0.,0.,0,LISO,0,6]
        self.head(tab6)
        tab6=[AWI,2.e8,1,0,10,2014]
        self.head(tab6)
        tab6=[0.,0.,0,0,NWD,NXC]
        self.head(tab6)

        for line in infofile:
            self.string(line)

        for line in content:
            tab6=['blnk','blnk',line[0],line[1],line[2],1]
            self.head(tab6)

        nline=self.IS+1
        self.send()
        self.fill_MF(451,nline)


class MF3(endf_section):
    """
    class description
    """

    def __init__(self,ZA=None,AWR=None,MAT=None):
        """
        function description
        """
        endf_section.__init__(self,ZA,AWR,MAT,3,None,0)
        self.filedata=[]
        self.ZA=ZA
        self.AWR=AWR
        self.MAT=MAT
        self.MF=3
        self.IS=1
        self.entry={}
        self.pos=0

    def set_MF(self,tab,entry):
        """
        function description
        """
        self.filedata=tab[:]
        self.entry=entry

    def MF_end(self):
        """
        function description
        """
        fendline= f' 0.000000+0 0.000000+0          0          0          0          0{self.MAT:4d} 0  0    0\n'
        self.filedata.append(fendline)
        MF_end= f'MF{self.MF:d}_end'
        self.entry[MF_end]=[self.pos,1]
        self.pos=self.pos+1

    def sort(self):
        """
        function description
        """
        MT_list=[]
        for key in self.entry.keys():
            i3=key.index('_')+3
            if key[i3] == 'd':
                MF3_end=self.entry['MF3_end']
                end_data=self.filedata[MF3_end[0]]
            else:
                i3=int(i3)
                MT=int(key[i3:])
                MT_list.append(MT)
        MT_list.sort()

        #store in temporary files
        entry_list=[]
        data=[]
        for i,MT in enumerate(MT_list):
            MF_MT='MF3_MT%i' %MT
            iMT=self.entry[MF_MT][0]
            nline=self.entry[MF_MT][1]

            entry_list.append([MF_MT,[i,nline]])
            data.append(self.filedata[iMT])

        #set sorted entry
        self.entry={}
        for t in entry_list:
            self.entry[t[0]] =t[1]
        self.entry['MF3_end']=MF3_end

        #set sorted data
        self.filedata=[]
        self.filedata.extend(data)
        self.filedata.append(end_data)


    def build_MT(self,MT,data):
        """
        function description
        """
        self.data=[]
        self.IS=0
        self.MT=MT
        QM=data[0]
        QI=data[1]
        LR=0
        tab6=[self.ZA,self.AWR,0,0,0,0]
        self.head(tab6)

        tab6=[QM,QI,0,LR,0,0]
        nbt=data[2][0]
        xs=data[2][1]
        self.tab1(tab6,nbt,xs)
        self.send()
        self.fill_MF(MT,self.IS)


    def get_QM(self,MT):
        """
        function description
        """
        MF_MT= f'MF3_MT{MT:d}'
        if MF_MT in self.entry:
            iMT=self.entry[MF_MT][0]
            line=self.filedata[iMT][1]
            QM=self.get_tab6(line)[0]
            return get_value(QM)

        print (f'section {MT:d} not present in file')
        return None

    def get_QValue(self,MT):
        """
        function description
        """
        MF_MT= f'MF3_MT{MT:d}'
        if MF_MT in self.entry:
            iMT=self.entry[MF_MT][0]
            line=self.filedata[iMT][1]
            QV=self.get_tab6(line)[1]
            return get_value(QV)

        print (f'section {MT:d} not present in file')
        return None

    def get_XS(self,MT,check=False):
        """
        function description
        """
        MF_MT= f'MF3_MT{MT:d}'
        if self.entry.has_key(MF_MT):
            iMT=self.entry[MF_MT][0]
            self.data=self.filedata[iMT]
            xs=self.get_tab1(1)
            if check:
                x0=0
                for i,xy in enumerate(xs[1][1]):
                    if (xy[0] < x0 ):
                        print (f' In MT{MT:d},  X value at {i:d} is not in increasing order')
                    x0=xy[0]
            return xs[1]
        else:
            print (f'section {MT:d} not present in file')
            return None

class MF6(endf_section):
    """
    class description
    """

    def __init__(self,ZA=None,AWR=None,MAT=None):
        """
        function description
        """
        endf_section.__init__(self,ZA,AWR,MAT,6,None,1)
        self.filedata=[]
        self.ZA=ZA
        self.AWR=AWR
        self.MAT=MAT
        self.MF=6
        self.IS=1
        self.entry={}
        self.pos=0

    def set_MF(self,tab,entry):
        """
        function description
        """
        self.filedata=tab[:]
        self.entry=entry

    def MF_end(self):
        """
        function description
        """
        fendline= f' 0.000000+0 0.000000+0          0          0          0          0{self.MAT:4d} 0  0    0\n'
        self.filedata.append(fendline)
        MF_end= f'MF{self.MF:d}_end'
        self.entry[MF_end]=[self.pos,1]
        self.pos=self.pos+1

    def seek_particule(self,MT,ZAP):
        """
        function description
        """
        MF_MT= f'MF6_MT{MT:d}'
        iMT=self.entry[MF_MT][0]
        line=self.filedata[iMT][0]
        tab6=self.get_tab6(line)
        npart=get_value(tab6[4])

        self.pointer=1
        part_init=1
        count_part=1

        Cont=True
        while Cont:
            if (count_part > npart ) :
                self.pointer=1
                print (f'Outgoing particle not found in MT{MT:d}')
                break
            line=self.filedata[iMT][self.pointer]
            tab6=self.get_tab6(line)
            part=get_value(tab6[0])
            nr  =get_value(tab6[4])
            np  =get_value(tab6[5])
            law =get_value(tab6[3])
            if (part==ZAP) :
                Cont=False
                part_init=self.pointer
            self.jump(2*nr,cont=True)
            self.jump(2*np)
            count_part=count_part+1
            if law == 1 or law == 2:
                line=self.filedata[iMT][self.pointer]
                tab6=self.get_tab6(line)
                nr  =get_value(tab6[4])
                nz  =get_value(tab6[5])
                self.jump(2*nr)
                for i in range(nz):
                    line=self.filedata[iMT][self.pointer]
                    tab6=self.get_tab6(line)
                    nl  =get_value(tab6[4])
                    self.jump(nl)
            elif law==7:
                line=self.filedata[iMT][self.pointer]
                tab6=self.get_tab6(line)
                nr  =get_value(tab6[4])
                n1  =get_value(tab6[5])
                self.jump(2*nr)
                for i in range(n1):
                    line=self.filedata[iMT][self.pointer]
                    tab6=self.get_tab6(line)
                    nr  =get_value(tab6[4])
                    n2  =get_value(tab6[5])
                    self.jump(2*nr)
                    for j in range(n2):
                        line=self.filedata[iMT][self.pointer]
                        tab6=self.get_tab6(line)
                        nr  =get_value(tab6[4])
                        n3  =get_value(tab6[5])
                        self.jump(2*nr,cont=True)
                        self.jump(2*n3)
            else:
                if ( law !=0 and law != 3 and law != 4):
                    print ('Unknown LAW number: ',law)


        self.part_len=self.pointer-part_init
        self.pointer=part_init
        self.part_point=self.pointer

    def remove_particle(self,MT,ZAP):
        """
        function description
        """

        self.seek_particule(MT,ZAP)
        if self.part_len == 0 :
            return

        self.MT=MT
        MF_MT= f'MF6_MT{MT:d}'
        iMT=self.entry[MF_MT][0]

        line= self.filedata[iMT][0]
        tab6=self.get_tab6(line)
        npart =get_value(tab6[4])
        tab6[4]=npart-1
        line = self.write_line(tab6)
        self.filedata[iMT][0]=line

        del self.filedata[iMT][self.part_point:self.part_point+self.part_len]
        pos=self.entry[MF_MT][0]
        nline=self.entry[MF_MT][1]-self.part_len
        self.entry[MF_MT]=[pos,nline]

    def copy_particle(self,MT,ZAP):
        """
        function description
        """

        self.seek_particule(MT,ZAP)
        if self.part_len == 0 :
            return

        self.MT=MT
        MF_MT= f'MF6_MT{MT:d}'
        iMT=self.entry[MF_MT][0]

        return self.filedata[iMT][self.part_point:self.part_point+self.part_len]

    def append_particle(self,MT,data):
        """
        function description
        """

        self.MT=MT
        MF_MT= f'MF6_MT{MT:d}'
        iMT=self.entry[MF_MT][0]

        line= self.filedata[iMT][0]
        tab6=self.get_tab6(line)
        npart =get_value(tab6[4])
        tab6[4]=npart+1
        line = self.write_line(tab6)
        self.filedata[iMT][0]=line

        self.filedata[iMT][-1:-1]=data
        pos=self.entry[MF_MT][0]
        nline=self.entry[MF_MT][1]+len(data)
        self.entry[MF_MT]=[pos,nline]
##########

    def get_particle(self,MT):
        """
        function description
        """
        MF_MT= f'MF6_MT{MT:d}'
        iMT=self.entry[MF_MT][0]
        line=self.filedata[iMT][0]
        tab6=self.get_tab6(line)
        npart=get_value(tab6[4])

        self.pointer=1
        #part_init=1
        count_part=1
        part_tab=[]

        Cont=True
        while Cont:
            if (count_part > npart ) :
                break
            line=self.filedata[iMT][self.pointer]
            tab6=self.get_tab6(line)
            part=get_value(tab6[0])
            nr  =get_value(tab6[4])
            np  =get_value(tab6[5])
            law =get_value(tab6[3])
            part_tab.append([part,law,self.pointer])

            self.jump(2*nr,cont=True)
            self.jump(2*np)
            count_part=count_part+1
            if law == 1 or law == 2:
                line=self.filedata[iMT][self.pointer]
                tab6=self.get_tab6(line)
                nr  =get_value(tab6[4])
                nz  =get_value(tab6[5])
                self.jump(2*nr)
                for i in range(nz):
                    line=self.filedata[iMT][self.pointer]
                    tab6=self.get_tab6(line)
                    nl  =get_value(tab6[4])
                    self.jump(nl)
            elif law==7:
                line=self.filedata[iMT][self.pointer]
                tab6=self.get_tab6(line)
                nr  =get_value(tab6[4])
                n1  =get_value(tab6[5])
                self.jump(2*nr)
                for i in range(n1):
                    line=self.filedata[iMT][self.pointer]
                    tab6=self.get_tab6(line)
                    nr  =get_value(tab6[4])
                    n2  =get_value(tab6[5])
                    self.jump(2*nr)
                    for j in range(n2):
                        line=self.filedata[iMT][self.pointer]
                        tab6=self.get_tab6(line)
                        nr  =get_value(tab6[4])
                        n3  =get_value(tab6[5])
                        self.jump(2*nr,cont=True)
                        self.jump(2*n3)
            else:
                if ( law !=0 and law != 3 and law != 4):
                    print ('Unknown LAW number: ',law)


        return part_tab

class MF8(endf_section):
    """
    class description
    """

    def __init__(self,ZA=None,AWR=None,MAT=None):
        """
        function description
        """
        endf_section.__init__(self,ZA,AWR,MAT,8,None,1)
        self.filedata=[]
        self.ZA=ZA
        self.AWR=AWR
        self.MAT=MAT
        self.MF=8
        self.IS=1
        self.entry={}
        self.pos=0

    def set_MF(self,tab,entry):
        """
        function description
        """
        self.filedata=tab[:]
        self.entry=entry

    def MF_end(self):
        """
        function description
        """
        fendline= f' 0.000000+0 0.000000+0          0          0          0          0{self.MAT:4d} 0  0    0\n'
        self.filedata.append(fendline)
        MF_end= f'MF{self.MF:d}_end'
        self.entry[MF_end]=[self.pos,1]
        self.pos=self.pos+1

    def get_halflife(self):
        """
        function description
        """

        if 'MF8_MT457' in self.entry:
            iMT=self.entry['MF8_MT457'][0]
            self.data=self.filedata[iMT]
            line=self.data[0]
            NST=self.get_tab6(line)[4]
            if (NST == 1) :
                return 'stable'
            else:
                line=self.data[1]
                T12=self.get_tab6(line)[0]
            return float(T12)
        else:
            print ('No MT457 section ')
            return None

    def get_gamma_spectra(self,Erg_min=0,extra=True):
        """
        function description
        """
        decay=[]


        if 'MF8_MT457' in self.entry:
            iMT=self.entry['MF8_MT457'][0]
            self.data=self.filedata[iMT]

            nrad_typ=self.get_tab6(self.data[0])[5]
            nc=self.get_tab6(self.data[1])[4]
            if nc == 6 :
                cline=3
            else:
                cline=8

            ndk=self.get_tab6(self.data[cline])[5]

            cline=int(cline+ndk)


            for it in range(nrad_typ):
                cline=int(cline+1)
                tab6=self.get_tab6(self.data[cline])
                rtyp=int(tab6[1])
                ner = tab6[5]

                if rtyp == 0 :
                    cline=cline+1
                    FD=self.get_tab6(self.data[cline])[0]
                    for ie in range(ner):
                        cline=int(cline+1)
                        tab6=self.get_tab6(self.data[cline])
                        Eg = tab6[0]
                        nt = int(tab6[4])
                        nloop=(nt-1)/6
                        cline=cline+1
                        Ri=self.get_tab6(self.data[cline])[2]
                        cline=cline+nloop
#                        if (Eg >= Emin and FD*Ri >= Imin):
                        if Eg*FD*Ri >= Erg_min :
                            decay.append([Eg,FD*Ri,0])

                elif (rtyp == 9 and extra ):
                    cline=cline+1
                    FD=self.get_tab6(self.data[cline])[0]
                    for ie in range(ner):
                        cline=cline+1
                        Eg=self.get_tab6(self.data[cline])[0]
                        cline=cline+1
                        Ri=self.get_tab6(self.data[cline])[2]
#                        if (Eg >= Emin and FD*Ri >= Imin):
                        if (Eg*FD*Ri >= Erg_min ):
                            decay.append([Eg,FD*Ri,9])

                else:
                    cline=cline+1+2*ner
        else:
            print ('No MT457 section ')

        return decay

    def get_avrE(self,styp='default',decay=True,detail=False):
        """
        function description
        """
        if decay:
            Edecay=self.get_avrE_decay()
            Elp=['Elp heat',Edecay[0]]
            Eem=['Eem heat',Edecay[1]]
            Ehp=['Ehp heat',Edecay[2]]
            Eavr=[Elp,Eem,Ehp]
            if (len(Edecay) > 3 and detail):
                name=['Eb-','Eb+','EAuger','ECIe','Egamma','EX-ray','EInB',\
                      'Eannih','Ealpha','Erecoil','Esf','En','Ep','Enu']
                for i in range(14):
                    Eavr.append([name[i],Edecay[i+3]])
        else:
            if styp == 'default' :
                styp=[0,9,1,2,4,6,8]
            elif styp == 'all':
                styp=[0,1,2,3,4,5,6,7,8,9]
            if not detail:
                EMlist=[]
                CPlist=[]
                for st in styp:
                    if st in [0,9]:
                        EMlist.append(st)
                    else:
                        CPlist.append(st)
                # MIO
                #EM_energy=get_EMlist_avrE(EMlist) # Creo que esto es un selt.get(abajo)
                #CP_energy=get_CPlist_avrE(CPlist)
                EM_energy = self.get_EM_avrE(EMlist)
                CP_energy = self.get_CP_avrE(CPlist)


                Eavr=[['Eem',EM_energy],['ECP',CP_energy]]
            else:
                Elist=[]
                Eavr=[]
                name=['gamma','b-','EC,b+','alpha','n','sf','p','e-','X-Ray']
                for st in styp:
                    if st in [0,9]:
                        energy=get_EM_avrE([st])
                    else:
                        energy=get_CP_avrE([st])
                    Elist.append([st,energy])

                for E in Elist:
                    Eavr.append([name[E[0]],E[1]])

        return Eavr

    def get_CP_avrE(self,slist=[1,2,4,7,8]):
        """
        function description
        """

        if 0 in slist:
            slist.remove(0)
        if 'MF8_MT457' in self.entry:
            iMT=self.entry['MF8_MT457'][0]
            self.data=self.filedata[iMT]

            nrad_typ=self.get_tab6(self.data[0])[5]
            nc=self.get_tab6(self.data[1])[4]
            if nc == 6 :
                cline=3
            else:
                cline=8

            ndk=self.get_tab6(self.data[cline])[5]
            cline=cline+ndk

            Eavr=0.
            for it in range(nrad_typ):
                cline=cline+1
                tab6=self.get_tab6(self.data[cline])
                styp=int(tab6[1])
                lcon=int(tab6[2])
                ner = tab6[5]

                # 1: beta rays; 2: b+ and/or EC; 4:alpha; 7: protons; 8: discrete e-
                if styp in slist :
                    cline=cline+1
                    Eavr = Eavr + self.get_tab6(self.data[cline])[2]
                    cline=cline+2*ner
                elif styp == 0 :
                    cline=cline+1
                    for _ in range(ner):
                        cline=cline+1
                        tab6=self.get_tab6(self.data[cline])
                        nt = int(tab6[4])
                        nloop=(nt-1)/6
                        cline=cline+nloop+1
                else:
                    cline=cline+1+2*ner
                if lcon != 0:
                    cline=cline+1
                    tab6=self.get_tab6(self.data[cline])
                    np=tab6[5]
                    nr=tab6[4]
                    lcov=tab6[3]
                    npl= (2*np-1)/6+1
                    nrl= (2*nr-1)/6+1
                    cline=cline+nrl+npl
                    if lcov != 0 :
                        cline=cline+1
                        tab6=self.get_tab6(self.data[cline])
                        npp=tab6[5]
                        nloop=(2*npp-1)/6
                        cline=cline+nloop+1

        else:
            print ('No MT457 section ')

        return Eavr


    def get_EM_avrE(self,slist=[0,9]):
        """
        function description
        """

        if 'MF8_MT457' in self.entry:
            iMT=self.entry['MF8_MT457'][0]
            self.data=self.filedata[iMT]

            nrad_typ=self.get_tab6(self.data[0])[5]
            nc=self.get_tab6(self.data[1])[4]
            if nc == 6 :
                cline=3
            else:
                cline=8

            ndk=self.get_tab6(self.data[cline])[5]
            cline=cline+ndk

            Eavr=0.
            for it in range(nrad_typ):
                cline=cline+1
                tab6=self.get_tab6(self.data[cline])
                styp=int(tab6[1])
                lcon=int(tab6[2])
                ner = tab6[5]

                if (styp == 0 ) :
                    cline=cline+1
                    if styp in slist:
                        Eavr = Eavr + self.get_tab6(self.data[cline])[2]
                    for ie in range(ner):
                        cline=cline+1
                        tab6=self.get_tab6(self.data[cline])
                        nt = int(tab6[4])
                        nloop=(nt-1)/6
                        cline=cline+nloop+1
                elif (styp == 9) :
                    cline=cline+1
                    if styp in slist:
                        Eavr = Eavr + self.get_tab6(self.data[cline])[2]
                    cline=cline+2*ner
                else:
                    cline=cline+1+2*ner
                if lcon != 0:
                    cline=cline+1
                    tab6=self.get_tab6(self.data[cline])
                    np=tab6[5]
                    nr=tab6[4]
                    lcov=tab6[3]
                    npl= (2*np-1)/6+1
                    nrl= (2*nr-1)/6+1
                    cline=cline+nrl+npl
                    if lcov != 0 :
                        cline=cline+1
                        tab6=self.get_tab6(self.data[cline])
                        npp=tab6[5]
                        nloop=(2*npp-1)/6
                        cline=cline+nloop+1

        else:
            print ('No MT457 section ')

        return Eavr

    def get_avrE_decay(self):
        """
        function description
        """

        if 'MF8_MT457' in self.entry:
            iMT=self.entry['MF8_MT457'][0]
            self.data=self.filedata[iMT]

            nrad_typ=self.get_tab6(self.data[0])[5]
            nc=self.get_tab6(self.data[1])[4]
            cline=2
            if nc == 6 :
                nl=1
            else:
                nl=6
            Elist=[]
            for i in range(nl):
                tab=self.get_tab6(self.data[cline])
                Elist.extend([tab[0],tab[2],tab[4]])
                cline=cline+1

            return Elist

            ndk=self.get_tab6(self.data[cline])[5]
            cline=cline+ndk

            for it in range(nrad_typ):
                cline=cline+1
                tab6=self.get_tab6(self.data[cline])
                styp=int(tab6[1])
                lcon=int(tab6[2])
                ner = tab6[5]

                if styp == 0 :
                    cline=cline+1
                    for ie in range(ner):
                        cline=cline+1
                        tab6=self.get_tab6(self.data[cline])
                        nt = int(tab6[4])
                        nloop=(nt-1)/6
                        cline=cline+nloop+1
                else:
                    cline=cline+1+2*ner
                if lcon != 0:
                    cline=cline+1
                    tab6=self.get_tab6(self.data[cline])
                    np=tab6[5]
                    nr=tab6[4]
                    lcov=tab6[3]
                    npl= (2*np-1)/6+1
                    nrl= (2*nr-1)/6+1
                    cline=cline+nrl+npl
                    if lcov != 0 :
                        cline=cline+1
                        tab6=self.get_tab6(self.data[cline])
                        npp=tab6[5]
                        nloop=(2*npp-1)/6
                        cline=cline+nloop+1

        else:
            print ('No MT457 section ')

        return Eavr


    def get_neutron_decay(self):
        """
        function description
        """

        if 'MF8_MT457' in self.entry:
            iMT=self.entry['MF8_MT457'][0]
            self.data=self.filedata[iMT]

            nsp=self.get_tab6(self.data[0])[5]
            nc=self.get_tab6(self.data[1])[4]
            if nc == 6 :
                cline=3
            else:
                cline=8

            ndk=self.get_tab6(self.data[cline])[5]
            for it in range(ndk):
                cline=cline+1
                rtyp=self.get_tab6(self.data[cline])[0]
                if int((rtyp+0.001)*10) in [15,50] :
                    br=self.get_tab6(self.data[cline])[4]
                    return [rtyp,br]

            for it in range(nsp):
                cline=cline+1
                tab6=self.get_tab6(self.data[cline])
                styp=int(tab6[1])
                lcon=int(tab6[2])
                ner = tab6[5]

                if (styp == 0) :
                    cline=cline+1
                    for ie in range(ner):
                        cline=cline+1
                        tab6=self.get_tab6(self.data[cline])
                        nt = int(tab6[4])
                        nloop=(nt-1)/6
                        cline=cline+nloop+1
                else:
                    cline=cline+1+2*ner
                if lcon != 0:
                    cline=cline+1
                    tab6=self.get_tab6(self.data[cline])
                    np=tab6[5]
                    nr=tab6[4]
                    lcov=tab6[3]
                    npl= (2*np-1)/6+1
                    nrl= (2*nr-1)/6+1
                    cline=cline+nrl+npl
                    if lcov != 0 :
                        cline=cline+1
                        tab6=self.get_tab6(self.data[cline])
                        npp=tab6[5]
                        nloop=(2*npp-1)/6
                        cline=cline+nloop+1

        else:
            print ('No MT457 section ')

        return None

    def get_alpha_decay(self):
        """
        function description
        """

        if 'MF8_MT457' in self.entry:
            iMT=self.entry['MF8_MT457'][0]
            self.data=self.filedata[iMT]

            nsp=self.get_tab6(self.data[0])[5]
            nc=self.get_tab6(self.data[1])[4]
            if nc == 6 :
                cline=3
            else:
                cline=8

            ndk=self.get_tab6(self.data[cline])[5]
            for it in range(ndk):
                cline=cline+1
                rtyp=self.get_tab6(self.data[cline])[0]
                if ( int((rtyp+0.001)*10) in [40,14,24] ) :
                    br=self.get_tab6(self.data[cline])[4]
                    return [rtyp,br]

            for it in range(nsp):
                cline=cline+1
                tab6=self.get_tab6(self.data[cline])
                styp=int(tab6[1])
                lcon=int(tab6[2])
                ner = tab6[5]

                if (styp == 0) :
                    cline=cline+1
                    for ie in range(ner):
                        cline=cline+1
                        tab6=self.get_tab6(self.data[cline])
                        nt = int(tab6[4])
                        nloop=(nt-1)/6
                        cline=cline+nloop+1
                else:
                    cline=cline+1+2*ner
                if lcon != 0:
                    cline=cline+1
                    tab6=self.get_tab6(self.data[cline])
                    np=tab6[5]
                    nr=tab6[4]
                    lcov=tab6[3]
                    npl= (2*np-1)/6+1
                    nrl= (2*nr-1)/6+1
                    cline=cline+nrl+npl
                    if lcov != 0 :
                        cline=cline+1
                        tab6=self.get_tab6(self.data[cline])
                        npp=tab6[5]
                        nloop=(2*npp-1)/6
                        cline=cline+nloop+1

        else:
            print ('No MT457 section ')

        return None


class MF12_13(endf_section):
    """
    class description
    """

    def __init__(self,ZA,AWR,MAT,MF):
        """
        function description
        """
        endf_section.__init__(self,ZA,AWR,MAT,MF,None,1)
        self.filedata=[]
        self.ZA=ZA
        self.AWR=AWR
        self.MAT=MAT
        self.MF=MF
        self.IS=1
        self.entry={}
        self.pos=0

    def set_MF(self,tab,entry):
        """
        function description
        """
        self.filedata=tab[:]
        self.entry=entry

    def MF_end(self):
        """
        function description
        """
        fendline= f' 0.000000+0 0.000000+0          0          0          0          0{self.MAT:4d} 0  0    0\n'
        self.filedata.append(fendline)
        MF_end='MF{self.MF:d}_end'
        self.entry[MF_end]=[self.pos,1]
        self.pos=self.pos+1

    def build_MT(self,MT,xstab):
        """
        function description
        """

        IS=self.IS
        self.MT=MT
        self.data=[]
        LO=1
        nk=len(xstab)
        tab6=[self.ZA,self.AWR,LO,0,nk,0]
        self.head(tab6)
        if nk == 1 :
            tab6 = xstab[0][0]
            nbt  = xstab[0][1]
            xs1  = xstab[0][2]
            self.tab1(tab6,nbt,xs1)
        else:
            tab6 = [0.,0.,0,0,0,0]
            nbt  = xstab[0][1]
            xs=[]
            for i in range(nk):
                xs.append(xstab[i][2])
            xstot=sum_tab1(xs)
            self.tab1(tab6,nbt,xstot)

            for xs in xstab:
                tab6 = xs[0]
                nbt  = xs[1]
                xs1  = xs[2]
                self.tab1(tab6,nbt,xs1)

        nline=self.IS-IS
        self.send()

        self.fill_MF(MT,nline)


class MF14(endf_section):
    """
    class description
    """

    def __init__(self,ZA,AWR,MAT):
        """
        function description
        """
        endf_section.__init__(self,ZA,AWR,MAT,14,None,1)
        self.filedata=[]
        self.ZA=ZA
        self.AWR=AWR
        self.MAT=MAT
        self.MF=14
        self.IS=1
        self.entry={}
        self.pos=0

    def set_MF(self,tab,entry):
        """
        function description
        """
        self.filedata=tab[:]
        self.entry=entry

    def MF_end(self):
        """
        function description
        """
        fendline= f' 0.000000+0 0.000000+0          0          0          0          0{self.MAT:4d} 0  0    0\n'
        self.filedata.append(fendline)
        MF_end= f'MF{self.MF:d}_end'
        self.entry[MF_end]=[self.pos,1]
        self.pos=self.pos+1

    def append_iso(self,MT,nk):
        """
        function description
        """

        IS=self.IS
        self.MT=MT
        self.data=[]

        LI=1   # isotropic for all photons

        tab6=[self.ZA,self.AWR,LI,0,nk,0]
        self.head(tab6)

        nline=self.IS-IS
        self.send()
        self.fill_MF(MT,nline)

class MF15(endf_section):
    """
    function description
    """

    def __init__(self,ZA,AWR,MAT):
        """
        function description
        """
        endf_section.__init__(self,ZA,AWR,MAT,15,None,1)
        self.filedata=[]
        self.ZA=ZA
        self.AWR=AWR
        self.MAT=MAT
        self.MF=15
        self.IS=1
        self.entry={}
        self.pos=0

    def set_MF(self,tab,entry):
        """
        function description
        """
        self.filedata=tab[:]
        self.entry=entry

    def MF_end(self):
        """
        function description
        """
        fendline= f' 0.000000+0 0.000000+0          0          0          0          0{self.MAT:4d} 0  0    0\n'
        self.filedata.append(fendline)
        MF_end= f'MF{self.MF:d}_end'
        self.entry[MF_end]=[self.pos,1]
        self.pos=self.pos+1

    def append(self,MT,cont_data):
        """
        function description
        """

        IS=self.IS
        self.MT=MT
        self.data=[]

        nc=len(cont_data)

        tab6=[self.ZA,self.AWR,0,0,nc,0]
        self.head(tab6)
        for subsec in cont_data:
            subhead=subsec[0]
            subdata=subsec[1]

            tab6=subhead[0]
            nbt=subhead[1][0]
            tabhead=subhead[1][1]
            self.tab1(tab6,nbt,tabhead)

            ne=len(subdata)
            tab6=[0,0,0,0,0,ne]
            nbt=subhead[2]
            self.tab2(tab6,nbt)
            for i,data in enumerate(subdata):
                nbt=data[0]
                tab=data[1]
                E=tabhead[i][0]
                tab6=[0,E,0,0,0,0]
                self.tab1(tab6,nbt,tab)

        nline=self.IS-IS
        self.send()
        self.fill_MF(MT,nline)

class ENDF:
    """
    class description
    """

    def __init__(self):
        """
        function description
        """
        self.file=[]
        self.entry={}
        self.file_nuc_index={}
        self.pos=0
        self.head_EAFMT=[]
        self.EAFMT_prnt=False
        self.ZAMAT=None

    def clear(self,data='data'):
        """
        function description
        """
        cleartab=['data']
        if data == 'all' :
            cleartab.append('index')

        if 'data' in cleartab:
            self.file=[]
            self.entry={}
            self.pos=0
            self.head_EAFMT=[]
            self.EAFMT_prnt=False

        if 'index' in cleartab:
            self.file_nuc_index={}

    def get_ZAMAT(self):
        """
        function description
        """
        if self.ZAMAT is not None :
            return self.ZAMAT[0],self.ZAMAT[1],self.ZAMAT[2]

    def open_endf(self,filename,EAFMT3=False):
        """
        function description
        """
        self.fic=open(filename)
        self.EAFMT3=EAFMT3

    def close_endf(self):
        """
        function description
        """
        self.fic.close()

    def mend(self):
        """
        function description
        """
        eofline='                                                                     0 0  0    0\n'
        self.file.append(eofline)

    def fend(self):
        """
        function description
        """
        eofline='                                                                    -1 0  0    0'
        self.file.append(eofline)

    def identify_nuclide(self,nuc,liso=0):
        """
        function description
        """
        if type(nuc) is str:
            i1=nuc.index('-')+1
            elemt=nuc[i1:i1+2].strip()
            i1=i1+3
            Am=int(nuc[i1:i1+3])
            i1=i1+3
            if nuc[i1:].strip()== '' :
                meta=''
            else:
                meta='m'
            nuclide= f'{elemt:s}{Am:d}{meta:s}'
        else:
            M=['','m','s','t']
            Z=int(nuc/1000)
            A=int(nuc-1000*Z)
            nuclide= f'{z_sym(Z):s}{A:d}{M[liso]:s}'
        return nuclide

    def generate_file_index(self):
        """
        function description
        """
        M=['','m','s','t']
        #read header
        line_len=81
        if self.EAFMT3 :
            line_len=76
        pos_init=self.fic.tell()
        cont= True
        while cont:
            pos=self.fic.tell()
            line=self.fic.readline()
            if len(line) >= line_len:
                if int(line[70:72]) != 0 :
                    self.fic.seek(pos,0)
                    cont=False

        #njump=5
        njump=1
        if (self.EAFMT3) :
            njump=3
        nextfile= True
        while nextfile:
            pos=self.fic.tell()
            for i in range(njump):
                line=self.fic.readline()

            if line == '':
                break
            current_nuc=line[0:11].strip()
            current_mat=int(line[66:70])
            if (current_mat == -1) :
                nextfile=False
                break

            if (not self.EAFMT3) :
                line=self.fic.readline()
            iso=int(line[34:44])
            if (iso > 3 ) :
                iso = 0

            #if (self.EAFMT3) : current_nuc=get_value(current_nuc)
            current_nuc=get_value(current_nuc)
            nuclide=self.identify_nuclide(current_nuc)

            while True:
                line=self.fic.readline()
                if line == '':
                    nextfile=False
                    break
                mat=int(line[66:70])
                if (mat == 0) :
                    break
            self.file_nuc_index[nuclide.lower()+M[iso]]=[pos,current_mat]
        self.fic.seek(pos_init,0)

    def get_list(self):
        """
        function description
        """
        return self.file_nuc_index.keys()

    def write_endf(self,filename):
        """
        function description
        """
        fic=open(filename,'w')

        for MF_file in self.file :
            for MT_sect in MF_file :
                for line in MT_sect :
                    fic.write(line)
        fic.close()

        if self.EAFMT_prnt:
            fname=filename+'.eafhead'
            self.write_head_EAFMT(fname)


    def read_endf(self,nuclide=None,rewind=True,EAFMT_prnt=False):
        """
        function description
        """
     #print('Read data',nuclide,rewind,EAFMT_prnt)

        if self.EAFMT3 :
            self.EAFMT_prnt=EAFMT_prnt
        self.file=[]
        self.entry={}
        iMF=0
        iMT=0
        nline=0
        MT_section=[]
        MF_file=[]
        if rewind:
            self.fic.seek(0,0)

        line_len=81
        if self.EAFMT3 :
            line_len=76

        #read header
        cont= True
        while cont:
            pos=self.fic.tell()
            line=self.fic.readline()
            if len(line) < line_len:
                MT_section.append(line)
            elif (line[66:70] in ['    ','   0','   1','  99']):
                MT_section.append(line)
                pos=self.fic.tell()
                cont=False
            else:
                self.fic.seek(pos,0)
                cont=False

        #jump to specified nuclide file
        #if no nuclide defined, take first of the list

        #njump=5
        njump=1
        if self.EAFMT3 :
            njump=3
        cont=True
        if nuclide is not None:
            if nuclide.lower() in self.file_nuc_index:
                position=self.file_nuc_index[nuclide.lower()][0]
                self.fic.seek(position,0)
            else:
                nextfile= True
                while nextfile:
                    pos=self.fic.tell()
                    for i in range(njump):
                        line=self.fic.readline()
                    nuc=line[0:11].strip()
                    line=self.fic.readline()
                    liso=int(line[33:44].strip())
                    nuc=get_value(nuc)
                    current_nuc=self.identify_nuclide(nuc,liso)
                    if current_nuc.lower() == nuclide.lower():
                        nextfile=False
                        self.fic.seek(pos,0)
                    else:
                        while True:
                            line=self.fic.readline()
                            if line == '':
                                print ('EOF reached. Nuclide not found')
                                nextfile=False
                                cont=False
                                break
                            mat=int(line[66:70])
                            if (mat == 0) :
                                break

        # read specified/first nuclide file
        pos_init=self.fic.tell()
        if self.EAFMT3 :
            self.fic.readline()
            self.fic.readline()
        line=self.fic.readline()

        eof=[]
        ZA=line[0:11]
        AWR=line[11:22]
        MAT=int(line[66:70])
        self.ZAMAT=[ZA,AWR,MAT]

        MT=0

        self.fic.seek(pos_init)

        skip2=False
        if self.EAFMT3  :
            skip2=True

        while cont :
            if skip2 :
                head1=self.fic.readline()
                head2=self.fic.readline()
                EAFMTlines=head1+head2
            line=self.fic.readline()
            if len(line)== 0 :
                cont=False
                break
            MAT=int(line[66:70])
            MF=int(line[70:72])
            MT=int(line[72:75])

            if MAT == 0:
                eof.append(line)
                tline='                                                                    -1 0  0    0'
                eof.append(tline)
                cont=False
                break

            if MF == 0 :
                MT_section.append(line)
                MF_file.append(MT_section)
                self.file.append(MF_file)
                MT_section=[]
                MF_file=[]

                MF_MT=MF_str+'end'
                self.entry[MF_MT]=[iMF,iMT,nline]
                iMF=iMF+1
                iMT=0
                nline=0
            elif MT == 0 :
                MT_section.append(line)
                MF_file.append(MT_section)
                MT_section=[]
                MF_MT = MF_str + MT_str
                self.entry[MF_MT]=[iMF,iMT,nline]
                iMT=iMT+1
                nline=0
                if self.EAFMT3  :
                    pos=self.fic.tell()
                    line=self.fic.readline()
                    self.fic.seek(pos)
                    MF=int(line[70:72])
                    if self.EAFMT_prnt:
                        self.head_EAFMT.append([ZA,MAT,mth,EAFMTlines])
                    if MF == 0:
                        skip2=False
                    else:
                        skip2=True
            else:
                MT_section.append(line)
                nline=nline+1
                MF_str='MF%i_' %MF
                MT_str='MT%i' %MT
                skip2=False
                mth=MT
        self.file.append(eof)

    def write_head_EAFMT(self,filename):
        """
        function description
        """
        fic=open(filename,'w')

        for reac in self.head_EAFMT:
            line=f'{reac[0]:s} {reac[1]:d} {reac[2]:d}\n'
            line=line+reac[3]
            fic.write(line)
        fic.close()


    def remove_MF(self,MF):
        """
        function description
        """

        MF_MT= f'MF{MF:d}_end'
        if MF_MT in self.entry:
            iMF=self.entry[MF_MT][0]
            del self.file[iMF]

            for key in self.entry.keys() :
                if (self.entry[key][0] == iMF ):
                    del (self.entry[key])
                elif (self.entry[key][0] > iMF ):
                    self.entry[key][0]=self.entry[key][0]-1

    def remove_MT(self,MF,MT):
        """
        function description
        """

        MF_MT=f'MF{MF:d}_MT{MT:d}'
        if MF_MT in self.entry:

            iMF=self.entry[MF_MT][0]
            iMT=self.entry[MF_MT][1]
            del (self.file[iMF][iMT])
            for key in self.entry.keys() :
                if (self.entry[key][0] == iMF ):
                    if (self.entry[key][1] == iMT ):
                        del (self.entry[key])
                    elif (self.entry[key][1] > iMT ):
                        self.entry[key][1]=self.entry[key][1]-1


            MF_end= f'MF{MF:d}_end'
            if self.entry[MF_end][1] == 0 :
                del (self.entry[MF_end])
                del (self.file[iMF][0])

    def copy_MF(self,MF):
        """
        function description
        """

        MF_MT= f'MF{MF:d}_end'
        if MF_MT in self.entry:
            iMF=self.entry[MF_MT][0]
            MT_entry={}
            for key in self.entry.keys():
                if self.entry[key][0]== iMF :
                    MT_entry[key]=[self.entry[key][1],self.entry[key][2]]
            return self.file[iMF],MT_entry
        else:
            return None,None


    def copy_MT(self,MF,MT):
        """
        function description
        """

        MF_MT= f'MF{MF:d}_MT{MT:d}'
        if MF_MT in self.entry:

            iMF=self.entry[MF_MT][0]
            iMT=self.entry[MF_MT][1]
            return self.file[iMF][iMT]
        else:
            return []


    def insert_MF(self,MF_file):
        """
        function description
        """

        MF_end= f'MF{MF_file.MF:d}_end'
        if MF_end in self.entry :
            print (f'MF file {MF_file.MF:d} already in the file')
            return

        MF_min=0
        iFM=-1
        for key in self.entry.keys():
            if key[3] == '_' :
                k=3
            else:
                k=4
            MF=int(key[2:k])

            if  MF_min < MF < MF_file.MF  :
                iFM=self.entry[key][0]
                MF_min=MF

        iFM=iFM+1
        self.file[iFM:iFM]=[MF_file.filedata]

        for key in self.entry.keys():
            if (self.entry[key][0] >= iFM ):
                self.entry[key][0] = self.entry[key][0]+1

        for key in MF_file.entry.keys():
            self.entry[key]=[iFM,MF_file.entry[key][0],MF_file.entry[key][1]]

    def get_Emax(self):
        """
        function description
        """
        if not 'MF1_MT451' in self.entry:
            print ('No MF information file to update')
            return None

        infofile=self.file[0][0]
        Emax=get_value(infofile[3][11:22])
        if (Emax < 1) :
            print (f'Emax in info file is {Emax:e}. Emax is set to 20 MeV')
            Emax=2e7
        return Emax

    def update_header(self,ZA,AWR,MAT):
        """
        function description
        """
        if not 'MF1_MT451' in self.entry:
            print ('No MF information file to update')
            return

        keys=sort_MFMT(self.entry)
        infofile=self.file[0][0]

        njump=int(infofile[4][44:55])
        nsect=len(keys)

        upinfo=infofile[0:njump+5]
        line4=upinfo[4][0:55]+endf_format(nsect)+upinfo[4][66:]
        upinfo[4]= line4

        self.entry['MF1_MT451'][2]=njump+nsect+4
        info=endf_section(ZA,AWR,MAT,1,451,njump+5)

        i=0
        for val in keys:
            key=val[0]
            MF =val[1]
            MT =val[2]
            line=infofile[njump+5+i]
            old=info.get_tab6(line)
            if (MF == old[2] and MT == old[3]):
                i=i+1
                mod=old[5]
            else:
                mod=0
            nline=self.entry[key][2]
            tab6=['blnk','blnk',MF,MT,nline,mod]
            upinfo.append(info.write_line(tab6))

        sendline= f' 0.000000+0 0.000000+0          0          0          0          0{MAT:4d} 1  099999\n'
        if upinfo[0][66:70] == '    ':
            firstline=upinfo[0][0:66]+'   0'+upinfo[0][70:]
            upinfo[0]=firstline
        upinfo.append(sendline)

        self.file[0][0]=upinfo

    def get_MF(self):
        """
        function description
        """
        MF_tab=[]
        for key in self.entry.keys():
            MF=int(key[2:key.index('_')])
            if MF not in MF_tab:
                MF_tab.append(MF)
        MF_tab.sort()
        return MF_tab

    def get_MT_old(self,MF):
        """
        function description
        """
        MT_tab=[]
        MF_str= f'MF{MF:s}_'
        for key in self.entry.keys():
            if MF_str in key:
                MT_str=key[key.index('_'):]
                if 'end' not in MT_str:
                    MT_tab.append(int(MT_str[3:]))
        MT_tab.sort()
        return MT_tab

    def get_MT(self,MF,MTsort=True):
        """
        function description
        """
        MT_tab=[]
        MF_str= f'MF{MF:s}_'

        if MTsort:
            for key in self.entry.keys():
                if MF_str in key:
                    MT_str=key[key.index('_'):]
                    if 'end' not in MT_str:
                        MT_tab.append(int(MT_str[3:]))
            MT_tab.sort()

        else:
            invert={}
            for MF_MT in self.entry.keys():
                if MF_str in MF_MT:
                    MT_str=MF_MT[MF_MT.index('_'):]
                    if 'end' not in MT_str:
                        invert[self.entry[MF_MT][1]]=int(MT_str[3:])
            tab=invert.keys()
            tab.sort()
            for key in tab:
                MT_tab.append(invert[key])

        return MT_tab

    def section_len(self,MF,MT):
        """
        function description
        """

        MF_MT= f'MF{MF:d}_MT{MT:d}'
        if MF_MT in self.entry:
            iMF=self.entry[MF_MT][0]
            iMT=self.entry[MF_MT][1]
            return len(self.file[iMF][iMT])
        else:
            return 0

#######################
