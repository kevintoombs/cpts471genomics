COMMANDS INPUT:
PS C:\Users\guenw\Documents\gitrepos\cpts471genomics\program1> .\Program1-Toombs\x64\Release\Program1-Toombs.exe .\Program1-Toombs\Program1-Toombs\ref\Human-Mouse-BRCA2-cds.fasta 1 .\Program1-Toombs\Program1-Toombs\ref\demo.config > newlocalout.txt
PS C:\Users\guenw\Documents\gitrepos\cpts471genomics\program1> .\Program1-Toombs\x64\Release\Program1-Toombs.exe .\Program1-Toombs\Program1-Toombs\ref\Human-Mouse-BRCA2-cds.fasta 0 .\Program1-Toombs\Program1-Toombs\ref\demo.config > newglobalout.txt

Fell free to run:
$Toombs-Program1.2>Program1-Toombs.exe Human-Mouse-BRCA2-cds.fasta 0 demo.config > out.txt
			The 0 for global alignment, or the 1 for local.
$Toombs-Program1.2>Program1-Toombs.exe Human-Mouse-BRCA2-cds.fasta 1 demo.config > out.txt

Output is found in newlocalout.txt and newglobalout.txt 
The previous ones are left just to show close, but slightly off, the previous implementation was.