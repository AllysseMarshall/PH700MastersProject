#
# upraveno
# 19.5. 2016, D. Trunec
#
##################################################################
def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""

##################################################################
def find_between_r( s, first, last ):
    try:
        start = s.rindex( first ) + len( first )
        end = s.rindex( last, start )
        return s[start:end]
    except ValueError:
        return ""

##################################################################
def form_subroutine(str_react):
  string="\
subroutine derivs(y,dydx)\n\
 IMPLICIT NONE\n\
 integer,  parameter :: DP=kind(1.0d0)\n\
 REAL(DP), DIMENSION(:), INTENT(IN) :: y\n\
 REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx\n\
"
  for i in range(len(str_react)):
    string+= str_react[i] + "\n"
  string += "end subroutine derivs\n"
  return string


##################################################################
def form_subroutine_radau(str_react):
  string="\
subroutine derivs(N,t,y,dydx,RPAR,IPAR)\n\
 IMPLICIT NONE\n\
 integer,  parameter :: DP=kind(1.0d0)\n\
 integer, INTENT(IN) :: N \n\
 REAL(DP),               INTENT(IN) :: t\n\
 REAL(DP), DIMENSION(N), INTENT(IN) :: y\n\
 REAL(DP), DIMENSION(N), INTENT(OUT) :: dydx\n\
 REAL(DP), DIMENSION(:), INTENT(IN) :: RPAR\n\
 INTEGER,  DIMENSION(:), INTENT(IN) :: IPAR\n\
"
  for i in range(len(str_react)):
    string+= str_react[i] + "\n"
  string += "end subroutine derivs\n"
  return string
##################################################################
def form_subroutine_vode(str_react):
  string="\
subroutine derivs(NEQ,T,y,dydx)\n\
 IMPLICIT NONE\n\
 integer,  parameter :: DP=kind(1.0d0)\n\
 INTEGER, INTENT (IN) :: NEQ\n\
 REAL(DP), INTENT (IN) :: T\n\
 REAL(DP), DIMENSION(NEQ), INTENT(IN) :: y\n\
 REAL(DP), DIMENSION(NEQ), INTENT(OUT) :: dydx\n\
"
  for i in range(len(str_react)):
    string+= str_react[i] + "\n"
  string += "end subroutine derivs\n"
  return string



##################################################################
def form_subroutine_JACOBI(string_jac,INTERFACE):

  if INTERFACE == "VODE":
   string="\n\n\n\n\
SUBROUTINE F_JAC(NEQ,T,Y,ML,MU,PD,NRPD)\n\
  IMPLICIT NONE\n\
!  INTEGER, INTENT (IN) :: NEQ, ML, MU, NRPD\n\
!  DOUBLE PRECISION, INTENT (IN) :: T\n\
!  DOUBLE PRECISION, INTENT (IN) :: Y(NEQ)\n\
!  DOUBLE PRECISION, INTENT (OUT) :: PD(NRPD,NEQ)\n\
  INTEGER :: NEQ, ML, MU, NRPD\n\
  DOUBLE PRECISION :: T\n\
  DOUBLE PRECISION :: Y(NEQ)\n\
  DOUBLE PRECISION :: PD(NRPD,NEQ)\n\
\n\
! PD(i,j) = \n\
\n\
"
  elif INTERFACE == "RADAU5":
   string="\n\n\n\n\
SUBROUTINE F_JAC(NEQ,T,Y,PD,NRPD,RPAR,IPAR)\n\
  IMPLICIT NONE\n\
  INTEGER, INTENT (IN) :: NEQ, NRPD\n\
  DOUBLE PRECISION, INTENT (IN) :: T\n\
  DOUBLE PRECISION, INTENT (IN) :: Y(NEQ)\n\
  DOUBLE PRECISION, INTENT (IN) :: RPAR(:)\n\
  INTEGER, INTENT (IN) :: IPAR(:)\n\
  DOUBLE PRECISION, INTENT (OUT) :: PD(NRPD,NEQ)\n\
\n\
! PD(i,j) = \n\
\n\
"
  else:
   string="\n\n\n\n\
SUBROUTINE F_JAC(NEQ,T,Y,ML,MU,PD,NRPD)\n\
  IMPLICIT NONE\n\
  INTEGER, INTENT (IN) :: NEQ, ML, MU, NRPD\n\
  DOUBLE PRECISION, INTENT (IN) :: T\n\
  DOUBLE PRECISION, INTENT (IN) :: Y(NEQ)\n\
  DOUBLE PRECISION, INTENT (OUT) :: PD(NRPD,NEQ)\n\
\n\
! PD(i,j) = \n\
\n\
"

#  for i in range(len(str_react)):
#    string+= str_react[i] + "\n"

#  for i in range(len(str_react)):

#     for j in range(len(str_react)):  
  for i in range(NumSpecies):
     for j in range(NumSpecies):        
        if string_jac[i][j] != "":
          string += "\n\n PD(%d,%d) = & \n" % (j+1,i+1) 
          string +=  string_jac[i][j] 
        else:
          string += "\n\n PD(%d,%d) = 0.0_DP" % (j+1,i+1)
  
  string += "\n RETURN \n\n" + "END SUBROUTINE F_JAC\n"
  return string


##################################################################

def form_rates( str_rates ):
  string="\
subroutine compute_rates( IPAR, RPAR, Rates_out )\n\
 IMPLICIT NONE\n\
! integer,  parameter :: DP=kind(1.0d0)\n\
 INTEGER , DIMENSION(:), INTENT(IN) :: IPAR\n\
 REAL(DP), DIMENSION(:), INTENT(IN) :: RPAR\n\
 REAL(DP), DIMENSION(1:NUM_REACTIONS), INTENT(OUT), optional  :: Rates_out\n\
"

  string +=  "\n" +  local_lines_rates + "\n"


  for i in range(len(str_rates)):
    string+= str_rates[i] + "\n"

#  string += "print *, rate\n"
  string += " if(present(Rates_out)) Rates_out = Rate \n"
  string += "end subroutine compute_rates\n"
  return string

##################################################################

###################################################################
#def expand_for_loop( line_tmp ):
#    expanded_string = ""
#    loop_str = find_between(line_tmp,"[","]")
#    line_to_expand = line_tmp.split("]")[1]
#    variable =  loop_str.split("=")[0]
#    iterace  =  loop_str.split("=")[1]
#    it_fin =  int(iterace.split(":")[1])
#    it_deb =  int(iterace.split(":")[0])
#    for it in range(it_deb,it_fin+1):
#       expanded=line_to_expand.replace("$"+variable, str(it))
#       subs_list = re.findall("(?<=\${)(.*?)(?=})",expanded)
#       for m in range(len(subs_list)):
#          expanded=expanded.replace("${"+subs_list[m]+"}",str(eval(subs_list[m])),1)
# 
#       expanded_string += expanded + "\n"
#    return expanded_string
# 
#######################################################################
#
#def exand_ForFor_loop( line_tmp ):
#    expanded_string = ""
#    line_to_expand = line_tmp.split("]")[1]    
#    loop_strings = find_between(line_tmp,"[","]")
##    print( loop_strings, loop_str_a, loop_str_b
##    loop_str_a = loop_strings.split(",")[0]
#
#    loop_str = loop_strings.split(",")
#    print( loop_str, len(loop_str)
#    sys.exit()
#   
#    variable_a =  loop_str_a.split("=")[0]
#    iterace_a  =  loop_str_a.split("=")[1]
#    it_fin_a =  int(iterace_a.split(":")[1])
#    it_deb_a =  int(iterace_a.split(":")[0]) 
#    print( variable_a, it_fin_a, it_deb_a
###    sys.exit()
##
#    loop_str_b = loop_strings.split(",")[1]
#    variable_b =  loop_str_b.split("=")[0]
#    iterace_b  =  loop_str_b.split("=")[1]
#    it_fin_b =  int(iterace_b.split(":")[1])
#    it_deb_b =  int(iterace_b.split(":")[0]) 
#    print( variable_b, it_fin_b, it_deb_b
#    for it_a in range(it_deb_a,it_fin_a+1):
#      for it_b in range(it_deb_b,it_fin_b+1):
#        expanded=line_to_expand.replace("$"+variable_a, str(it_a)).replace("$"+variable_b, str(it_b))
#        subs_list = re.findall("(?<=\${)(.*?)(?=})",expanded)
#        for m in range(len(subs_list)):
#          expanded=expanded.replace("${"+subs_list[m]+"}",str(eval(subs_list[m])),1)
#
#        expanded_string += expanded + "\n"
#        print( expanded_string
#    return expanded_string
#


## TRICK to implement for for any number of iterations:
#
#  import itertools
#  pool=[{} for i in range(3)]
#  pool[0]=range(4);
#  pool[1]=range(3);
#  pool[2]=range(2);
#  print( list(itertools.product(*pool)) # list of all iterations
##

def expand_for_loop(line_tmp):
  expanded_string = ""
  line_to_expand = line_tmp.split("]")[1]
  loop_strings = find_between(line_tmp,"[","]")
  loop_str_parts =  loop_strings.split(",")
  num_loops = len(loop_str_parts)
  pool=[{} for i in range(0,num_loops)]
  variables = [ i.split("=")[0] for i in loop_str_parts]
  tmp_ranges    = [ i.split("=")[1] for i in loop_str_parts]  
  it_fin = [ i.split(":")[1] for i in tmp_ranges ]
  it_deb = [ i.split(":")[0] for i in tmp_ranges ]
  for it in range(0,num_loops):
    pool[it] = range(int(it_deb[it]),int(it_fin[it])+1)
  indexes_list = list(itertools.product(*pool))

  for it in range(len(indexes_list)):
    expanded=line_to_expand
    for var in range(num_loops):
#      print( iters, variables[var]
      expanded=expanded.replace("$"+variables[var], str(indexes_list[it][var]))
    subs_list = re.findall("(?<=\${)(.*?)(?=})",expanded)
#    print( subs_list # expanstion of expressions for examle ${$i+3} = 4 with $j=1
    for m in range(len(subs_list)):
      expanded=expanded.replace("${"+subs_list[m]+"}",str(eval(subs_list[m])),1)

    expanded_string += expanded + "\n"

  return expanded_string

########################################################################
#


import sys
import re
import itertools

f=open(sys.argv[1],'r')

# Remove comments form input file
data=""
for line in f:
  line = line.split("#",1)[0].strip()
  #expand for loops:
  if line.startswith("for"):
     data +=  expand_for_loop(line) + "\n"
#  elif line.startswith("ForFor"):
#     data += expand_ForFor_loop(line) + "\n"
#  elif line.startswith("FOR"):
#     data += expand_FOR_loop(line) + "\n"
  else:
     data+=line+"\n"


f=open("react_tmp","w")
f.write(data)
f.close()



#  data+=line.split("#",1)[0].strip()+"\n"

###############################
# get the interface INTERFACE = VODE ? RK4
interface_str =  find_between( data, "INTERFACE:","ENDINTERFACE")
interface_str = filter(None,interface_str.split("\n"))
INTERFACE = interface_str[0]
############################
  
reactions_str = find_between( data, "REACTIONS:", "ENDREACTIONS")

reactions_str = filter(None,reactions_str.split("\n"))


##################################
# add local variables to compute_rates subroutine
local_lines_rates = "" 
for line in reactions_str[:]:
  if line.startswith("$$"):
     local_lines_rates +=  line.split("$$",1)[1].strip() + "\n"
     reactions_str.remove(line)
###################################


##############################################################
species_str = find_between( data, "SPECIES:", "ENDSPECIES" )
fixed_species_str = find_between( data, "FIXED:", "ENDFIXED")


species_str =  filter(None,species_str.replace("\n"," ").strip().split(" "))


#reactions_str = filter(None,reactions_str.split("\n"))

#for line in reactions_str:
#  if line.startswith()

fixed_str = filter(None,fixed_species_str.replace("\n"," ").strip().split(" "))
indxs_fixed =  [species_str.index(item) for item in fixed_str]


##########################################################################
# bulid list of species
list_of_species=""
for line in species_str[:]:
   list_of_species +=  "! %-4d " % (species_str.index(line) + 1)
   list_of_species +=  line + "\n"



NumSpecies   = len(species_str)
NumReactions = len(reactions_str)


##############################################################
initvals_str =  find_between(data,"INITVAL:", "ENDINITVAL")

init_vals_subroutine = "\n\n\
subroutine initiate_reaction_scheme(DENSITIES, IPAR, RPAR )\n\
 IMPLICIT NONE \n\
 real(DP), intent(OUT), dimension(:) :: DENSITIES\n\n\
 integer, intent(IN), dimension(:)  :: IPAR\n\
 real(DP), intent(IN), dimension(:) :: RPAR\n\n\
\
"


initvals_str = filter(None, initvals_str.split("\n"))

### extract local code after "$$"
for line in initvals_str[:]:
  if line.startswith("$$"):
     init_vals_subroutine +=  line.split("$$",1)[1].strip() + "\n"
     initvals_str.remove(line)
##   print( line
##
for line in initvals_str[:]:
   split = filter(None,line.strip().split(" "))
   if split[0] == "DEFAULT": # DEFAULT
      init_vals_subroutine += "\n DENSITIES = " + "".join(split[1:]) + "   ! set DEFAULT values"  + "\n\n"
   else: # IDENTIFY species
      init_vals_subroutine += " DENSITIES(%d) = " % (species_str.index(split[0])+1) + "".join(split[1:]) + "\n"
  

# initiation of Jacobian
#init_vals_subroutine +=  "\n PD_JAC = 0.0_DP \n\n"



#init_vals_subroutine += DENSITIES 

#print( initvals_str
#print( init_vals_subroutine

#exit
############################################################


##################################################################

react_arr = [[0]*NumSpecies for i in range(NumReactions)]

RATE = [0]*NumReactions
RATE_VAL = [0]*NumReactions
indxs_reactants = [None]*NumReactions


for I_reaction in range(NumReactions):
   REAC = reactions_str[ I_reaction ]     
   REACTION =  REAC.split("!")[0] 

   rate_str_tmp = REAC.split("!")[1]

   print (REACTION)
#### BOLSIG
   if rate_str_tmp.split()[0] == "BOLSIG":
     list_tmp = rate_str_tmp.split()
     part_name_tmp = list_tmp[1].replace("(", " ").replace(")"," ").split()
     if len(list_tmp) == 3:
#        print( "bolsig_%s_%04d_( en )*%s" % (part_name_tmp[0],int(part_name_tmp[1] ,list_tmp[2]) )
        BOLSIG_STRING  = "bolsig_%s_%04d_( EN )*%s ! %s" % (part_name_tmp[0],int(eval(part_name_tmp[1])) ,list_tmp[2], "bolsig_en: e/N [Td]" )
        print (BOLSIG_STRING)
     else:
        BOLSIG_STRING = "bolsig_%s_%04d_( EN )" % (part_name_tmp[0],int(eval(part_name_tmp[1])))
        print (BOLSIG_STRING)

     RATE_VAL[I_reaction]     = BOLSIG_STRING

#### BOLSIG Sumire
   if rate_str_tmp.split()[0] == "bolsig":
     list_tmp = rate_str_tmp.split()
#     print( len(list_tmp), list_tmp
#     part_name_tmp = list_tmp[1].replace("(", " ").replace(")"," ").split()
     if len(list_tmp) == 4:
       BOLSIG_STRING = "bolsig_%s_%s( EN ) %s" % (list_tmp[1],list_tmp[2],list_tmp[3])       
##        print( "bolsig_%s_%04d_( en )*%s" % (part_name_tmp[0],int(part_name_tmp[1] ,list_tmp[2]) )
#        BOLSIG_STRING  = "bolsig_%s_%s_( EN )*%s ! %s" % (part_name_tmp[0],int(eval(part_name_tmp[1])) ,list_tmp[2], "bolsig_en: e/N [Td]" )
#        print( BOLSIG_STRING
     else:
       BOLSIG_STRING = "bolsig_%s_%s( EN )" % (list_tmp[1],list_tmp[2])
     print( BOLSIG_STRING)

     RATE_VAL[I_reaction]     = BOLSIG_STRING



#### POPOV
   elif rate_str_tmp.split()[0] == "POPOV":
     list_tmp = rate_str_tmp.split()
     if len(list_tmp) == 5:
#        print(  "POPOV" 
        params =  [ "%+.3e" % (float(x)) + "_dp" for x in list_tmp[1:-1]]
        popov_rate_str_tmp = "10**(%%s  %%s/%s  %%s/%s**2 )"  % (list_tmp[-1], list_tmp[-1] ) 
        popov_rate_str =  popov_rate_str_tmp % tuple(params)
        
     else:
        print(  list_tmp)
        print( " parsereac.py: parameter missing or too many\n\n")
        exit()
     
     RATE_VAL[I_reaction]     =  popov_rate_str

   elif rate_str_tmp.split()[0] == "POPOV_STAR":
     list_tmp = rate_str_tmp.split()
     if len(list_tmp) == 5:
#        print(  "POPOV_STAR")
        order_tmp = ["1e-9","1e-11","1e-13"]
#        print(  ["{}*{}".format(b_, a_) for a_, b_ in zip(order_tmp,list_tmp[1:-1])])
        params =  [ "%+.3e" % (eval(x)) + "_dp" for x in ["{}*{}".format(b_, a_) for a_, b_ in zip(order_tmp,list_tmp[1:-1])] ]
        popov_rate_str_tmp = "%%s  %%s*%s  %%s*%s**2"  % (list_tmp[-1], list_tmp[-1] ) 
        popov_rate_str =  popov_rate_str_tmp % tuple(params)
        
     else:
        print(  list_tmp)
        print( " parsereac.py: parameter missing or too many\n\n")
        exit()
     
     RATE_VAL[I_reaction]     =  popov_rate_str




   else:
     RATE_VAL[I_reaction]     =  REAC.split("!")[1]


   RATE[I_reaction] = "Rate(%d)" % (I_reaction+1)
   SIDES = REACTION.split("=>")
   
   reactants =  SIDES[0].strip().split(" + ")
   reactants = [item.strip() for item in reactants]

   products  =  SIDES[1].strip().split(" + ")
   products  = [item.strip() for item in products]
   
   indxs_reactants[I_reaction] = [species_str.index(item) for item in reactants]
   indxs_products  = [species_str.index(item) for item in products]

  
   NumReactants = len(indxs_reactants[I_reaction])
   NumProducts = len(indxs_products)
     
   for i in indxs_reactants[I_reaction]:
      react_arr[I_reaction][i] -= 1 
   
   for i in indxs_products:
      react_arr[I_reaction][i] += 1 
 
#   print(  reactions_str[ I_reaction ]
#   print(  [ x+1 for x in indxs_reactants[I_reaction]], "=> ",  [ x+1 for x in  indxs_products] ,"#", I_reaction
#   print( " $$$$$ ", react_arr[I_reaction][:], "\n"

       

 
#########
rate_set_val_str = [""]*NumReactions

for I_reaction in range(NumReactions):
  rate_set_val_str[I_reaction] = "Rate(%d) = " % (I_reaction+1) + RATE_VAL[I_reaction] 
#  print( rate_set_val_str[I_reaction]
#########################################################3

code_spec_line=[""]*NumSpecies
lhs = [ "\n dydx(%d) = & \n " % (i+1) for i in range(NumSpecies) ]


## jacobian
str_jac = [x[:] for x in [[""]*NumSpecies]*NumSpecies]


for j in range(NumReactions): 
      code_str = ""
      for k in indxs_reactants[j]:
         code_str += "y(" + str(k+1) + ")*"

      code_str += RATE[j]

      for i in range(NumSpecies):

#### indent long lines, if to be added to them
         if react_arr[j][i] != 0:
           if len(code_spec_line[i].split("&")[-1]) > 100 :
             code_spec_line[i] += " & \n"
      
         if react_arr[j][i] == 1:
              code_spec_line[i] +=  " + " + code_str

         elif react_arr[j][i] == -1:
              code_spec_line[i] +=  " - " + code_str

         elif react_arr[j][i] != 0:
            code_spec_line[i] +=    "%+d" % react_arr[j][i] + "*" + code_str


##############
# Jacobian
#      print(  reactions_str[ j ]
#      print(  [ x+1 for x in indxs_reactants[j]], "=> ",  [ x+1 for x in  indxs_products] ,"#", j
#      print( " $$$$$ ", react_arr[j][:], "\n"

      for i in range(NumSpecies): 
#         j: reaction, i specie
         if react_arr[j][i] != 0:
#              print(  [ x+1 for x in indxs_reactants[j]], i , "<--"
              for k in indxs_reactants[j]:
#                  print( "    ", indxs_reactants[j], k
                  tmp_list = list(indxs_reactants[j])
                  tmp_list.remove(k)
                  if react_arr[j][i] == 1:
                    tmp_str = " +"
                  elif react_arr[j][i] == -1:
                    tmp_str = " -"
                  else:
                    tmp_str = " %+d*" % react_arr[j][i]
                  for r in tmp_list:
                     tmp_str += "y(%d)*" % (r+1)                  
                  tmp_str += "Rate(%d)" % (j+1)
#                  print( tmp_str
                
                  if len(str_jac[k][i].split("&")[-1]) > 100 :
                     str_jac[k][i] += "&\n "

                  str_jac[k][i] += tmp_str
                  
#      print( str_jac
#                  tmp_list = indxs_reactants[j]
#                  for r in tmp_list:
#                     tmp_list = indxs_reactants[j]
#                     tmp_list.remove(r)
#                     print( "    y(%d)" % (r+1),  indxs_reactants[j]
#                 print(  " PD(%d,%d) =  " % (i+1,k+1)
#                 for r in indxs_reactants[j]: 
#                   tmp_list = indxs_reactants[j]
#                  tmp_list.remove(r)
#                   print(  "  & \n %+d" %  react_arr[j][i], indxs_reactants[j], r
#                   for l in tmp_list[:]:  
#                      print( "    y(%d)" % (l+1)
                 

# , "y(%d) = " % k ,   indxs_reactants[j], i , "*",react_arr[j][i],"*", k
#             for k in range(len(indxs_reactants[j])):
#                print( k,  
#      exit()


###### Those species not evolved:
for i in range(NumSpecies):
  if code_spec_line[i] == "":
     code_spec_line[i] = " 0.0_dp"

### Those are fixed by definition:
for i in indxs_fixed:
     code_spec_line[i] = " 0.0_dp"
#######################################

out = "\
module reaction_scheme \n "
out += find_between( data, "MODULES:", "ENDMODULES")
out +="\
use params\n\
implicit none\n\
!integer, parameter :: dp = kind(1.0d0)\n\
public :: compute_rates, derivs, initiate_reaction_scheme, NUM_SPECIES, F_JAC, NUM_REACTIONS\n\
private\n\
integer, parameter :: NUM_SPECIES = %d \n\
integer, parameter :: NUM_REACTIONS = %d \n\
REAL(DP), DIMENSION(1:NUM_REACTIONS)  ::  Rate\n\
!$OMP THREADPRIVATE(Rate) \n\
!REAL(DP), DIMENSION(1:NUM_SPECIES,1:NUM_SPECIES)  ::  PD_JAC\n\
" % (NumSpecies,NumReactions) \
 + find_between( data, "DECLARATIONS:", "ENDDECLARATIONS") + "\n contains \n"


for i in range(NumSpecies):
   code_spec_line[i] = lhs[i] + code_spec_line[i]

if INTERFACE == 'VODE':
  out += form_subroutine_vode(code_spec_line)
elif INTERFACE == 'RADAU5' :
  out += form_subroutine_radau(code_spec_line)
else:
  out += form_subroutine(code_spec_line)


out +=  form_subroutine_JACOBI(str_jac,INTERFACE)
out += "\n\n\n" +  form_rates( rate_set_val_str )
out += init_vals_subroutine + "end subroutine initiate_reaction_scheme\n\n"
out += "end module reaction_scheme\n\n"
out += list_of_species

#print( out

f=open("reaction_scheme.f90","w")
f.write(out)
f.close()

######################################

