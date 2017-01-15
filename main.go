package main

import (
	//"bufio"
	"flag"
	"fmt"
	//"io/ioutil"
	"log"
	"math/rand"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

var addr string = ""
var motifs [][]int
var mutations []int
var sequences []string

func main() {
	fmt.Println("Use flags to run the program with parameters.  If you did not provide flags, the program will run with default values.  The proper flags are documented in the readme.pdf file.")
	//flag values
	l := flag.Int("length", 500, "an int")
	q := flag.Int("quantity", 50, "an int")
	s := flag.String("sizes", "20, 20, 30", "a string")
	n := flag.String("mutations", "3, 3, 4", "a string")
	m := flag.Int("minimum", 1, "an int")
	b := flag.Int("subregion", 100, "an int")
	flag.Parse()
	//fmt.Println(*l)
	max, motifSizes := maximal(*s)
	_, mutationCount := maximal(*n)
	mutations = mutationCount
	//check the values
	valid := true
	if *l < *b {
		valid = false
	} else if ((*m) * (max)) > *b {
		valid = false
	} else if *l < 1 || *q < 1 || max < 0 || len(mutations) < 0 || *m < 0 || *b < 0 {
		valid = false
	} else if len(mutations) != len(motifSizes) {
		valid = false
	}
	if valid {
		mo := makeMotifs(motifSizes)
		motifs = mo
		seqs := make(chan string, *q)
		mots := make(chan string, len(motifSizes))
		generate(*l, *q, *m, *b, seqs)
		//fmt.Println("Should be writing")
		writer(seqs, mots, *q, len(motifSizes))
	} else {
		fmt.Println("Invalid parameters detected. Program ending.")
	}
}

//this function converts a string to ints, and finds the largest of these ints
func maximal(s string) (int, []int) {
	strings := strings.SplitN(s, ", ", -1)
	sizes := make([]int, len(strings), len(strings))
	max := 0
	for i := 0; i < len(strings); i++ {
		t, err := strconv.Atoi(strings[i])
		if err != nil {
			return -1, sizes
		}
		sizes[i] = t
		if sizes[i] > max {
			max = sizes[i]
		}
	}
	return max, sizes
}

//this function generates the random data and sequences using concurrency
func generate(leng, quan, mini, subr int, seqs chan<- string) {
	//fmt.Println("gen")
	//fmt.Println(motifs)
	gr := make(chan []int)
	gmot := make(chan []int)
	kill := make(chan int)
	go genMotifs(gmot, kill)
	go genRandom(gr, leng, kill)
	sequences = make([]string, quan)
	for i := 0; i < quan; i++ {
		sequences[i] = write(i, leng, mini, subr, gr, gmot)
	}
	kill <- 1
	close(gr)
	close(gmot)
	//fmt.Println("ending generate")
	/*
		for true {
			//stall so the sub goroutines won't go out of scope. Main will end while this function stalls, completing the tasks
		}
	*/

}

//this function generates the motifs, passes them into a channel
func genMotifs(gmots chan<- []int, kill <-chan int) {
	i := 0
	killed := false
	for !killed {
		select {
		case <-kill:
			{
				killed = true
				break
			}
		default:
			break
		}
		if killed {
			break
		}
		i = rand.Intn(len(motifs))
		t := make([]int, 1, 1)
		t[0] = i
		gmots <- append(t, mutate(i)...)
	}
}

//this function generates random noise and passes it into a channel
func genRandom(gr chan<- []int, size int, kill <-chan int) {
	s, i := size/10, 0
	killed := false
	for !killed {
		select {
		case <-kill:
			{
				killed = true
				break
			}
		default:
			break
		}
		if killed {
			break
		}
		i = rand.Intn(s)
		t := make([]int, i)
		for j := range t {
			t[j] = rand.Intn(4) + 1
		}
		gr <- t
	}
}

//creates the motifs
func makeMotifs(sizes []int) [][]int {
	motifs := make([][]int, len(sizes), len(sizes))
	for i := 0; i < len(motifs); i++ {
		temp := make([]int, sizes[i], sizes[i])
		for i := range temp {
			temp[i] = rand.Intn(4) + 1
			// 1-4 represent real nucleotides, 0 represents unfilled space later
		}
		motifs[i] = temp
		//fmt.Println(motifs[i])
	}
	//fmt.Println(motifs)
	return motifs
}

//writes the sequences
func write(ind, leng, mini, subr int, gr, gm <-chan []int) string {
	//fmt.Println("write called")
	str := "> seq" + strconv.Itoa(ind)
	//channels
	seq := make([]int, leng, leng)
	if len(motifs) == 0 {
		//fmt.Println("motifs list empty")
	} else {
		i := rand.Intn(leng - subr)
		j, k := i+subr, 0
		//fmt.Println("i ", i, " j ", j)

		//generate the subregion using channels
		for true {
			//fmt.Println("i ", i, " j ", j)
			i, k = j-subr, 0
			for i < j {
				//fmt.Println("motif len ", len(gm), " random len ", len(gr))
				select {
				case m := <-gm:
					{
						//fmt.Println("motif")
						i += len(m)
						if i > j {
							break
						}
						k++
						seq = append(seq[:i-len(m)], append(m[1:], seq[i:]...)...)
						str += " m" + strconv.Itoa(m[0]) + " " + strconv.Itoa(i-len(m)) + " " + strconv.Itoa(i)
					}
				case n := <-gr:
					{
						//fmt.Println("nonmotif")
						i += len(n)
						if i > j {
							break
						}
						//fmt.Println("i", i, "len(n)", len(n))
						seq = append(seq[:i-len(n)], append(n, seq[i:]...)...)
					}
				}
			}
			//fmt.Println("K: ", k)
			if k >= mini {
				break
			}
		}
	}
	//fill in the nonsubregion part of the sequence
	fill(&seq)
	str += "\n" + translate(seq) + "\n"
	//fmt.Println("writing: ", str)
	return str
	//fmt.Print(str)
	//if channel okays, you write to file/pass the seq through a channel to the writer
}

//single nucleotide mutations occur here
func mutate(i int) []int {
	muta := mutations[i]
	motif := motifs[i]
	taken := make([]int, len(motifs[i]))
	index, dice := 0, 0
	for i := 0; i < muta; i++ {
		index = rand.Intn(len(motif))
		if taken[index] == 1 {
			i--
			continue
		}
		taken[index] = 1
		dice = rand.Intn(3)
		if dice == 0 {
			//insertion
			motif = append(append(motif[:index+1], rand.Intn(4)+1), motif[index+1:]...)
			taken = append(append(taken[:index+1], 1), taken[index+1:]...)
		} else if dice == 1 {
			//deletion
			motif = append(motif[:index], motif[index+1:]...)
			taken = append(taken[:index], taken[index+1:]...)
		} else {
			motif[index] = rand.Intn(4) + 1
		}
	}
	return motif
}

//fill nonmutated regions randomly.
func fill(s *[]int) {
	seq := *s
	for i := range seq {
		if seq[i] == 0 {
			seq[i] = rand.Intn(4) + 1
		}
	}
	s = &seq
}

//unused function from older version
func pass(m chan<- string) {
	str := ""
	//fmt.Println("MOTIFS COUNT ", len(motifs))
	for i := range motifs {
		//fmt.Println("Passing mot ", i)
		str = translate(motifs[i])
		m <- str
	}
	close(m)
	//fmt.Println("closed mots")
}

//converts integers to nucleotides
func translate(seq []int) string {
	str := ""
	for _, i := range seq {
		switch i {
		case 1:
			str += "A"
		case 2:
			str += "C"
		case 3:
			str += "G"
		case 4:
			str += "T"
		}
	}
	return str
}

//writes to the output files
func writer(seqs, mots <-chan string, ms, ss int) {
	//fmt.Println("WRITER")
	dir, err := filepath.Abs(filepath.Dir(os.Args[0]))
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println("writing output files to ", dir)
	fs, err := os.Create("sequences.txt")
	if err != nil {
		fmt.Println("Error creating sequence file")
	}
	fm, err := os.Create("motifs.txt")
	if err != nil {
		fmt.Println("Error creating motif file")
	}
	//make two files
	//fmt.Print("mot", motif, mc, "\n")
	motif := ""
	for i, m := range motifs {
		motif = "m" + strconv.Itoa(i) + " " + translate(m) + "\n"
		_, err := fm.WriteString(motif)
		if err != nil {
			fmt.Println("Error writing motif file")
		}
	}
	for _, sequence := range sequences {
		_, err := fs.WriteString(sequence)
		if err != nil {
			fmt.Println("Error writing sequence file")
		}
	}
	fm.Close()
	fs.Close()
	fmt.Println("Finished writing")
}
