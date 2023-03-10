Ανάπτυξη Λογισμικού για Αλγοριθμικά Προβλήματα - Εργασία 3

Αλέξανδρος Θέμελης sdi1900062
Δημήτριος Γεωργαντόπουλος sdi1900036

Μεταγλώττιση: Όντας στο directory, καλούμε τις εντολές “cmake -DCGAL_DIR=usr/include/CGAL .” και ύστερα “make”. Έτσι, μπορεί πλέον να εκτελεστεί μέσω της υποδεδειγμένης εντολής της εκφώνησης.

ΣΗΜΑΝΤΙΚΟ:
Για εμάς, η προεπιλογή σημαίνει να δίνονται 3 σημειοσύνολα ανά αριθμό σημείων, και ορίζεται από τη σταθέρα "NUM_OF_FILES". Αν δοθούν λιγότερα ή περισσότερα σημειοσύνολα ανά αιρθμό σημείων, τότε να αλλαχθεί στη γραμμή 36 στο αρχείο functions.h.

Παραδείγματα εκτέλεσης:
./evaluate -i points -o output.txt

Λίστα αρχείων:
CMakeLists.txt (έχει επεξεργαστεί καταλλήλως για separate compilation)
evaluate_polygon.cpp (η main συνάρτηση, που διαβάζει, κάνει χρονομετρήσεις και βγάζει αποτελέσματα)
includes/CMakeLists.txt (έχει επεξεργαστεί καταλλήλως για separate compilation)
includes/convex_hull_alg.cpp (αλγόριθμος βάσει ΚΠ)
includes/functions.cpp (βοηθητικές συναρτήσεις για διάβασμα σημείων και area από input files,και sort)
includes/functions.h
includes/incremental_alg.cpp (αυξητικός αλγόριθμος)
includes/incremental_alg_subdivisial.cpp (αυξητικός αλγόριθμος παραλλαγμένος για χρήση του στο spatial subdivision)
includes/localSearch.cpp (Τοπική Αναζήτηση)
includes/simulatedAnnealing.cpp (Προσομοιωμένη ανόπτηση με υλοποιημένα τα local/global step και spatial subdivision)
readme.txt

Προτείνεται στο CMakeCache.txt που θα δημιουργηθεί να αλλαχθεί το CMAKE_BUILD_TYPE:STRING= στην γραμμή 57 ως CMAKE_BUILD_TYPE:STRING=Release για να τρέχουν πιο γρήγορα τα αποτελέσματα. 

Η έρευνα και δοκιμές έγιναν από κοινού με την ταυτόχρονη συνεργασία και των δυο μας σε έναν υπολογιστή, εξού και τα commit μόνο από ένα account. Και οι δύο μας έχουμε γνώση όλων των σημείων του προγράμματος και μπορούμε να εξηγήσουμε την εκτέλεσή του λεπτομερώς.

ΣΧΟΛΙΑΣΜΟΣ ΑΠΟΤΕΛΕΣΜΑΤΩΝ: 
1.incremental local search
Ο αυξητικός αλγόριθμος με τυχαία επιλογή ορατής ακμής χρησιμοποιείται στον αλγόριθμο καθώς είχε παρατηρηθεί από την πρώτη εργασία οτι παράγει αρκετά ικανοποιητικά αποτελέσματα με γρήγορη εκτέλεση. Την χρειαζόμαστε ιδιαίτερα, καθώς η τοπική αναζήτηση αποτελεί brute force αλγόριθμο και συνεπώς χρειάστηκε να ισοσταθμίσουμε πολυπλοκότητα και ταχύτητα. Φυσικά όμως, σε μεγάλα σημειοσύνολα, η brute force λογική προκαλεί την μεγάλη καθυστέρηση του αλγορίθμου. Γι'αυτό, όσο και να αλλάξουν τα ορίσματα του, δεν είναι αναμενόμενο ο αλγόριθμος να πληροί τις χρονικές προυποθέσεις για σημειοσύνολα άνω των 1000 σημείων. Δυναμικά αλλάζουμε το L και το threshold, με τιμές που έχουμε ερευνήσει ενδελεχώς, ώστε να καταφέρουμε όσο περισσότερη ακρίβεια και να είμαστε εντός ορίων του cutoff. Ο αλγόριθμος θεωρείται ικανοποιητικός, αλλά δεν μπορεί να είναι η βασική μας επιλογή, αφού μας έχουν ζητηθεί έλεγχος σε μεγάλο αριθμό σημείων.

2.incremental simulated annealing with local step
Εδώ πέρα χρησιμοποιήθηκε ο αυξητικός αλγόριθμος με simulated annealing με τοπικό βήμα, καθώς είχαμε παρατηρήσει ήδη από την δεύτερη εργασία, οτι ήταν από τους καλύτερους συνδυασμούς τεχνικών. Είναι γρήγορες και οι δυο στρατηγικές. Όσο μεγαλώνει το μέγεθος των σημειοσυνόλων, παρατηρείται αύξηση στον χρόνο εκτέλεσης, όχι τόσο ραγδαίο όμως όπως του incremental local search. Και σε αυτόν τον αλγόριθμο έγινε σημαντικός πειραματισμός με τις τιμές του L, καθώς προέκυψε οτι παίζει καθορισιτκό ρόλο στην αποτελεσματικότητα του αλγορίθμου.

3.convex hull/incremental simulated annealing with global step
Ο convex hull αλγόριθμος ενώ δείχνει σε μικρά σημειοσύνολα αρκετά καλά αποτελέσματα σε γρήγορο σχετικά χρόνο για την φύση του, σε μεγαλύτερα παρουσιάζει το μεγάλο του μειονέκτημα. H περιορισμένη επίδοσή του, ακόμα και με καθολικό βήμα, το πιο γρήγορο optimisation που έχουμε δηλαδή, δεν αφήνει να γενικευτούν τα αποτελέσματα που προέκυπταν πρίν. Συνεπώς, για σημειοσύνολα άνω των 2000 σημείων, ο αλγόριθμος θα ήταν εδώ πέρα για λόγους πληρότητας, καθώς οι υψηλές απαιτήσεις της άσκησης, δεν αφήνουν να φανούν οι καλές επιδόσεις του αλγορίθμου.Για τον λόγο αυτό, μετά από τα 1000 σημεία, η αρχική πολυγονοποίηση, αλλάζει και γίνεται από τον incremental. Έτσι, με υβριδικό, δύναται να επιτευχθούν οι θετικές πτυχές και των δυο αλγορίθμων που υλοποιήσαμε στη πρώτη εργασία. Εννοείται πως και εδώ χρησιμοποιήθηκε τυχαία η επιλογή ορατής ακμής, καθώς οταν δοκιμάστηκαν οι πιο εξειδικευμένες επιλογές, χειροτέρευσε κατά πού η σχέση ακρίβειας/αποδοτικότητας. Ο αλγόριθμος συνιστάται για σημειοσύνολα κάτω των 2000 σημείων. 

σημείωση:
Μπορεί το καθολικό βήμα να είναι ελαφρώς γρηγορότερο, αλλά συγκριτικά, η μεγαλύτερη ακρίβεια  αποτελεσμάτων που προσφέρει ο αλγόριθμος 2 αξίζει περισσότερο αν ασχοληθεί κανείς ιδιαιτέρως με τη τιμή του L. Συνολικά, προτιμάται περισσότερο ο αλγόριθμος 2.

Ο spatial subdivision ήταν υπερβολικά ασταθής ως αλγόριθμος, γι'αυτό και παραλείπεται ως επιλογή. Ακόμα, παρατηρήθηκε οτι τα αποτελέσματα που έδινε, δεν απέδιδε καλύτερα από τον βασικό μας αλγόριθμο, αλλά ήταν μέτριος, τόσο σε θέμα χρόνου, όσο και σε ακρίβεια. Ενώ στη θεωρία, η προσέγγισε διαίρει και βασίλευε ακούγεται ελκυστική, στη πράξη παρατηρήσαμε οτι σε μεγάλα σημειοσύνολα, αυξάνονται τα edge cases, δημιουργούνται δυσκολίες, και συνεπώς περιορίζεται η αποδοτικότητα του αλγορίθμου. Πιθανώς με μερικές υποθέσεις, ο αλγόριθμος θα τα πήγαινε καλύτερα.
