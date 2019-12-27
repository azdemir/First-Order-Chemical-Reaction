︠7cd1e86e-0c44-452b-b19c-50390969e5e8s︠
# Gas Absorption with First Order Reaction in a Packed Column with Surface-Renewal and Film Models. Full Models: both gas side and liquid side resistances are present. It is assumed that the concentration of dissolved gas in the inlet liquid is zero. The tower is assumed to be filled with 1 in Raschig rings. It is assumed that there is no dissovled gas in the liquid entering the tower at the top.
# Input parameters for air_SO2_H2O System at 21 C. Some of the parameters given below were calculated in a separate spreadsheet.
h = 5 # Packed height of tower (m)
d = 0.5 # Column diameter (m)
LM = 2.14 # Liquid mass velocity in tower (kg/m2/s)
RHOL = 998.204 # Liquid density at 20 C (kg/m3)
GM = 0.935 # Vapor mass velocity in tower (kg/m2/s)
RHOG = 1.2064 # Vapor density at 20 C (kg/m3)
H = 0.024304771369379 # Henry's constant at 21 C for the air-SO2-H2O system. Calculated in a separate spreadsheet.
HL_PA = 0.575106662952067 # Height of liquid-phase transfer unit for physical absorption (m). Calculated in a separate spreadsheet.
HG_PA = 1.0879503114122 # Height of gas-phase transfer unit for physical absorption (m) Calculated in a separate spreadsheet.
print "HL_PA (m) = ", HL_PA, " HG_PA (m) = ", HG_PA
kL = 0.000160344918705 # Liquid-phase mass-transfer coefficient (m/s). Estimated by the correlation of Onda et al. (1967) in a separate spreadsheet.
kG = 0.030642229077144 # Gas-phase mass-transfer coefficient (m/s). Estimated by the correlation of Onda et al. (1967) in a separate spreadsheet.
a = 23.2482820852841 # Gas-liquid interfacial area of contact (m2/m3). It is assumed that this is equal to the packing wetted area estimated by the correlation of Onda et al. (1967) in a separate spreadsheet.
phi = 0.036998817347497 # Fractional liquid hold-up (m3/m3). Calculated from eqs. (18-49), (18-51) and (18-52) in Perry's ChE Handbook (6th edition).
Diff = 1.49*10^(-9) # Diffusion coefficient of SO2 in H2O at 20 C (m2/s)
#Z = 0.215299641518421 # kL/(H kG)   Dimensionless quantity where H is Henry's Law coefficient and kL and kG are liquid- and gas-phase mass transfer coefficients (m/s)
Ha = 0.3 # Hatta number
print "Ha = ", Ha
# Start of Calculations
# SR Model
Z = kL/H/kG # Dimensionless parameter
L = LM/RHOL # Liquid superficial velocity in tower (m3/s/m2)
G = GM/RHOG # Vapor superficial velocity in tower (m3/s/m2)
E = sqrt(1 + Ha^2) # Enhancement factor
HOL = HL_PA/E + (L/H/G)* HG_PA # Height of overall liquid-phase transfer unit for chemical absorption (m)
HOG = HG_PA + (H*G/L)*HL_PA/E # Height of overall gas-phase transfer unit for chemical absorption (m)
print "HOG (SR, m) = ", HOG, " HOL (SR, m) = ", HOL
NOG = h/HOG # Number of gas-phase transfer units for chemical absorption
NOL = h/HOL # Number of liquid-phase transfer units for chemical absorption
print "NOG (SR, m) = ", NOG, " NOL (SR, m) = ", NOL
P = kL*phi/Diff/a  # Dimensionless parameter
K = P*Ha^2*(Z + 1/E) # Dimensionless reaction rate constant. Valid for both surface-renewal and film models.
f = 1/(1 + Ha^2)# Ratio of the rate of dissolved gas transfer to the bulk liquid to the rate of gas absorption at at the gas-liquid interface.
rsur(x) = f
Term1 = NOG - NOL*(K + f/E^2)
Term2 = sqrt((NOG + NOL*(K + f/E^2))^2 - 4*f*NOL*NOG/E^2)
lambda_A = 0.5*(Term1 + Term2)
lambda_B = 0.5*(Term1 - Term2)
yout = (lambda_B - lambda_A)/((NOG - lambda_A)*exp(lambda_B) - (NOG - lambda_B)*exp(lambda_A)) # Dimensionless outlet gas concentration at the top of the tower
print "yout (SR) = ", yout
DA = E^2*(1 - lambda_A/NOG)
DB = E^2*(1 - lambda_B/NOG)
c1 = -DB*yout/(DA - DB)
c2 = DA*yout/(DA - DB)
cout = c1*DA*exp(lambda_A) + c2*DB*exp(lambda_B) # Dimensionless outlet liquid concentration of dissolved gas at the bottom of the tower
print "cout (SR) = ", cout
# x is the dimensionless distance measured from the top of the column, x = 0 is the top of the column and x = 1 is the bottom of the column.
ysur (x) = c1*exp(lambda_A*x) + c2*exp(lambda_B*x) # Dimensionless vapor-phase concentration profile
csur(x) = c1*DA*exp(lambda_A*x) + c2*DB*exp(lambda_B*x) # Dimensionless liquid-phase  concentration profile
# Film Model
E = Ha/tanh(Ha) # Enhancement factor
HOL = HL_PA/E + (L/H/G)* HG_PA # Height of overall liquid-phase transfer unit for chemical absorption (m)
HOG = HG_PA + (H*G/L)*HL_PA/E # Height of overall gas-phase transfer unit for chemical absorption (m)
print "HOG (film, m) = ", HOG, " HOL (film, m) = ", HOL
NOG = h/HOG # Number of gas-phase transfer units for chemical absorption
NOL = h/HOL # Number of liquid-phase transfer units for chemical absorption
print "NOG (film, m) = ", NOG, " NOL (film, m) = ", NOL
f1 = 1/cosh(Ha)
E1 = sqrt(cosh(Ha))
E2 = 1/sqrt(cosh(Ha)*(1 + Z*Ha*tanh(Ha)))
Term1 = NOG - NOL*(K + f1/E2^2)
Term2 = sqrt((NOG + NOL*(K + f1/E2^2))^2 - 4*f1*NOL*NOG/E1^2)
lambda_A = 0.5*(Term1 + Term2)
lambda_B = 0.5*(Term1 - Term2)
yout = (lambda_B - lambda_A)/((NOG - lambda_A)*exp(lambda_B) - (NOG - lambda_B)*exp(lambda_A)) # Dimensionless outlet gas concentration at the top of the tower
print "yout (Film) = ", yout
DA = E1^2*(1 - lambda_A/NOG)
DB = E1^2*(1 - lambda_B/NOG)
c1 = -DB*yout/(DA - DB)
c2 = DA*yout/(DA - DB)
cout = c1*DA*exp(lambda_A) + c2*DB*exp(lambda_B) # Dimensionless outlet liquid concentration of dissolved gas at the bottom of the tower
print "cout (Film) = ", cout
# x is the distance measured from the top of the column, x = 0 is the top of the column and x = 1 is the bottom of the column.
yfilm(x) = c1*exp(lambda_A*x) + c2*exp(lambda_B*x) # Dimensionless vapor-phase concentration profile
cfilm(x) = c1*DA*exp(lambda_A*x) + c2*DB*exp(lambda_B*x) # Dimensionless liquid-phase concentration profile
rfilm(x) = f1*(yfilm(x) - H*cfilm(x)/E2^2)/(yfilm(x) - H*cfilm(x)/E1^2) # Ratio of the rate of gas transfer to the bulk liquid to the rate of gas absorption at the gas-liquid interface
# Create an empty plot object
B = plot([],figsize=(5,5),axes_labels=['dimensionless distance (z*)','dimensionless concentration (y*, c*)'],gridlines=True, frame=true, xmax = 1, ymax = 1)
B += plot(ysur(x),x,[0,1],color='blue', legend_label='SR Model (gas-phase)', linestyle="solid", thickness=1)
B += plot(csur(x),x,[0,1],color='red', legend_label='SR Model (liquid-phase)', linestyle="dashed", thickness=1)
B += plot(yfilm(x),x,[0,1],color='green', legend_label='Film Model (gas-phase)', linestyle="dashdot", thickness=1)
B += plot(cfilm(x),x,[0,1],color='black', legend_label='Film Model (liquid-phase)', linestyle="dotted", thickness=1)
B += text("Ha = 0.3 (LP and GP control)", (0.5,0.1), alpha=0.8, fontsize='large', fontweight='regular', color='black', axis_coords='True')
show(B)
#######
# Create an empty plot object
B = plot([],figsize=(5,5),axes_labels=['dimensionless distance (z*)','ratio of transfer to absorption rate (f)'],gridlines=True, frame=true, xmin=0, ymin=0, xmax = 1, ymax = 1)
B += plot(rsur(x),x,[0,1],color='blue', legend_label='SR Model', linestyle="solid", thickness=1)
B += plot(rfilm(x),x,[0,1],color='red', legend_label='Film Model', linestyle="dashed", thickness=1)
B += text("Ha = 0.3 (LP and GP control)", (0.5,0.1), alpha=0.8, fontsize='large', fontweight='regular', color='black', axis_coords='True')
show(B)
︡73841e09-bfc9-4924-8f06-23d16d29a6a0︡{"stdout":"HL_PA (m) =  0.575106662952067  HG_PA (m) =  1.08795031141220\n"}︡{"stdout":"Ha =  0.300000000000000\n"}︡{"stdout":"HOG (SR, m) =  5.92802992295193  HOL (SR, m) =  0.674672536949749\n"}︡{"stdout":"NOG (SR, m) =  0.843450533311444  NOL (SR, m) =  7.41100271044886\n"}︡{"stdout":"yout (SR) =  0.446455612751428\n"}︡{"stdout":"cout (SR) =  0.0482028120011089\n"}︡{"stdout":"HOG (film, m) =  5.99481166983960  HOL (film, m) =  0.682273006444715\n"}︡{"stdout":"NOG (film, m) =  0.834054558403464  NOL (film, m) =  7.32844470288325\n"}︡{"stdout":"yout (Film) =  0.451760726896677\n"}︡{"stdout":"cout (Film) =  0.0497997850061306\n"}︡{"file":{"filename":"/home/user/.sage/temp/project-78f81c40-5f91-441d-a7eb-ea999ed56ef8/655/tmp_2xkMx4.svg","show":true,"text":null,"uuid":"b670bcc9-fbc5-4594-80ab-fc03dc126e35"},"once":false}︡{"file":{"filename":"/home/user/.sage/temp/project-78f81c40-5f91-441d-a7eb-ea999ed56ef8/655/tmp_POMw7K.svg","show":true,"text":null,"uuid":"090ba56c-723f-4cbf-b3ef-eda9b3f098f6"},"once":false}︡{"done":true}









