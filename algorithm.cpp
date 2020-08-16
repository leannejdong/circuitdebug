#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <stack>
#include <iterator>
#include <cassert>
#include <fstream>
#include <sstream>


using namespace std;
using std::vector;

class Graph {
    int r;
    std::vector<std::vector<bool>> adjMatrix;
    std::vector<std::vector<int>> treeAdjMat;

    int countDifference(auto &treeAdjMat_i, auto &adjMatrix_i)
    {
        int n_differences = 0;
        for(size_t j=0; j<treeAdjMat_i.size(); j++)
        {
            if(treeAdjMat_i[j]!=adjMatrix_i[j])
                ++n_differences;
        }

        return n_differences;
    }

public:
    // Initialize the matrix to zero
    Graph(int r): r(r), adjMatrix(r, std::vector<bool>(r, false)),
                          treeAdjMat(r, std::vector<int>(r)) {}

    void addEdge(int i, int j) {
        assert(i >= 0 && i < r && j > 0 && j < r);
        adjMatrix[i][j] = true;
        adjMatrix[j][i] = true;
    }
    void removeEdge(int i, int j) {
        assert(i >= 0 && i < r && j > 0 && j < r);
        adjMatrix[i][j] = false;
        adjMatrix[j][i] = false;
    }
    bool isEdge(int i, int j) {
        if (i >= 0 && i < r && j > 0 && j < r)
        {
            return adjMatrix[i][j];
        }
        else
        {
            return false;
        }
    }

    const vector<vector<int>>&getTreeAdjMat() const
    {
        return treeAdjMat;
    }


    vector<vector<int>> Gotlieb(int r, int *m, vector<vector<bool>>& adjMatrix)
    {
        *m = 0; int i, j, k, c, nu, done, f, n;
        vector<int> x(r, 0);
        //Block 1
        vector<vector<int>> treeAdjMat(r, vector<int>(r));
        vector<vector<int>> copytreeAdjMat(r, vector<int>(r));
        for (i=0; i < r; i++)
        {
            done = 0;
            for (j=i; j < r; j++)
            {
                if(adjMatrix[i][j]==1 && done==0){
                    treeAdjMat[i][j]=1;
                    treeAdjMat[j][i]=1;
                    done=1;
                }
            }
        }
        //Block 2
        for(i=0; i<r; i++)
        {
            for(j=0; j<r; j++)
            {
                assert(i < r);
                copytreeAdjMat[i][j]=treeAdjMat[i][j];
            }
        }
        cout << r << "\n";

        for(i=0; i<r; i++)
        {
            for(j=0; j<r; j++)
            {
                if(copytreeAdjMat[i][j]==1)
                {
                    for(k=0; k<r; k++)
                        if(copytreeAdjMat[j][k]==1){
                            copytreeAdjMat[i][k]=copytreeAdjMat[j][k];
                            copytreeAdjMat[j][k]=-copytreeAdjMat[j][k];
                        }
                }
            }
        }

        // now join together the strands that are part of a single cord
        for(j=0; j<r; j++){
            bool found = false;
            for(i=0; i<r; i++)
            {
                if(copytreeAdjMat[i][j]==1 && !found)
                {
                    k=i;
                    found=true;
                }
                else if(copytreeAdjMat[i][j]==1 && found){
                    for(int n=0; n<r; n++)
                        if(copytreeAdjMat[i][n]==1)
                        {
                            copytreeAdjMat[k][n]=copytreeAdjMat[i][n];
                            copytreeAdjMat[i][n]=-copytreeAdjMat[i][n];
                        }
                }
            }
        }

        //count how many lines the matrix c has
        n=0;
        for(i=0;i<r; i++)
        {
            for(j=0;j<r; j++)
            {
                if(copytreeAdjMat[i][j]==1){
                    n++;
                    break;
                }
            }
        }

        vector<vector<int>> dmat(n, vector<int>(r));

        k=0;
        for(i=0; i<r; i++)
        {
            for(j=0; j<r; j++)
            {
                if(copytreeAdjMat[i][j]==1)
                {
                    for(f=0; f<r; f++)
                    {
                        dmat[k][f]=copytreeAdjMat[i][f];
                    }
                    k++;
                    break;
                }

            }
        }

        //correctness check
        for(j=0; j<r; j++)
        {
            for(i=0; i<n; i++)
            {
                assert(i<n);
                if(dmat[i][j]==1)
                {
                    continue;
                }
                if(dmat[i][j]==0 && j==n-1)
                {
                    printf("\nError in block 2 while searching for independent meshes\n");
                }
            }
        }

        // So far, we have constructed a spanning tree. Now let focus on finding cycle bases
        /*BLOCK 3*/
        f=0;
        for(i=0;i<n; i++)
            for(j=0;j<r; j++)
                if(dmat[i][j]==1)
                    for(k=0; k<r; k++) /*I'm going to see if I find a side that joins the cords in the j-th column of a*/
                        if(adjMatrix[k][j]==1 && dmat[i][k]==0 && f==0) {
                            /*printf("k=%d\tj=%d\n", k, j);*/
                            adjMatrix[k][j]=1;
                            adjMatrix[j][k]=1;
                            f=1;       /*I have to add only one side!*/
                        }

        for(i=0; i<r; i++)
        {
            vector<int> &treeAdjMat_i = treeAdjMat[i];
            vector<bool> &adjMatrix_i = adjMatrix[i];

            //int n_differences = 0;
            assert(static_cast<size_t>(r)==treeAdjMat_i.size());
            *m += countDifference(treeAdjMat_i,adjMatrix_i);
        }
        int &count = *m;
        count /= 2;
        //count how many sides have to be eliminated to obtain the tree graph = number of independent links
        c = r*count + count + 1;
        vector<vector<int>> indm(r);
        for (int i = 0; i<r; ++i)
        {
            indm[i].resize(c);
        }

        for (j = 0; j < c-r; j = j+r+1)
            for (i = 0; i < r; i++)
                indm[i][j] = -4;
        for (i = 0; i < r; i++)
            indm[i][c-1]=-5;
        for (k = 1; k < c; k=k+r+1)
            for(i = 0; i < r; i++)
                for(j = 0; j < r; j++)
                    indm[i][j+k] = treeAdjMat[i][j];
        // add the sides at a time
        k = 1;
        for(i = 0; i < r; i++)
            for(j = i+1; j<r; j++)
                if(adjMatrix[i][j]==1 && treeAdjMat[i][j]==0)
                {
                    indm[i][j+k]=1;
                    indm[j][i+k]=1;
                    k = k + r + 1;
                }
        /*I remove the segments that are outside the loop (see drawing)*/
        nu = 0; /*nu is the number one on a line*/
        done=0;
        for(k=1; k<c; k=k+r+1){
            while(done==0){
                done=1;
                for(i=0; i<r; i++){
                    for(j=0; j<r; j++)
                    {
                        /*Count how many ones are on a line*/
                        if(indm[i][j+k]==1)
                        {
                            nu++;
                        }
                    }
                    if(nu==1)
                    {
                        /*if there is only one,  make it null*/
                        for(j=0; j<r; j++)    /*I am in the j of 1*/
                        {
                            if(indm[i][j+k]==1){
                                indm[i][j+k]=0;
                                indm[j][i+k]=0;
                                done=0;
                            }
                        }
                    }
                    nu=0;
                }
            }
            done=0;
        }
        return indm;
    }

    void printMat() {
        int i, j;
        for (i = 0; i < r; i++ )
        {
            for (j = 0; j < r; j++)
            {
                std::cout << to_string(adjMatrix[i][j]) << " ";
            }
            std::cout << "\t";

            for (j = 0; j < r; j++)
            {
                std::cout << to_string(treeAdjMat[i][j]) << " ";
            }
            cout << endl;
        }
    }
};

// Requires a sequence of closed cycles.
template <class ForwardIterator, class OutputStream>
void print_cycles(ForwardIterator first, ForwardIterator last, OutputStream &os)
{
    using T = typename iterator_traits<ForwardIterator>::value_type;
    while (first != last)
    {
        auto const cycle_start = first++;
        first = ++find(first, last, *cycle_start);
        copy(cycle_start, first, ostream_iterator<T>(os, " "));
        os << "\n";
    }
}


static const char *circuit_matrix_text =
        "0 1 0 0 0 0 1 1 0 0 0\n"
        "1 0 1 0 0 0 0 0 0 0 0\n"
        "0 1 0 1 0 0 0 1 0 0 0\n"
        "0 0 1 0 1 0 0 1 0 0 0\n"
        "0 0 0 1 0 1 0 0 0 0 0\n"
        "0 0 0 0 1 0 1 0 0 0 0\n"
        "1 0 0 0 0 1 0 0 0 0 1\n"
        "1 0 1 0 0 0 0 0 1 0 0\n"
        "0 0 0 0 1 0 0 1 0 1 0\n"
        "0 0 0 0 0 0 0 0 1 0 1\n"
        "0 0 0 0 0 0 1 0 0 1 0\n";


static void printIndmTo(std::ostream &stream, const vector<vector<int>> &indm)
{
    for (size_t i=0; i!=indm.size(); ++i) {
        for (size_t j=0; j!=indm[i].size(); ++j) {
            if (j != 0) {
                stream << "  ";
            }
            stream << indm[i][j];
        }
        stream << "\n";
    }
}

static void printIndm(const vector<vector<int>> &indm)
{
    std::ostream &stream = std::cerr;

    for (size_t i=0; i!=indm.size(); ++i) {
        for (size_t j=0; j!=indm[i].size(); ++j) {
            if (j != 0) {
                stream << "  ";
            }
            stream << indm[i][j];
        }
        stream << "\n";
    }
}

static vector<vector<bool>> makeBoolMatrixFromIntMatrix(const vector<vector<int>> &indm)
{
    vector<vector<bool>> result(indm.size());
    for (size_t i=0; i!=indm.size(); ++i) {
        result[i].resize(indm[i].size());
        for (size_t j=0; j!=indm[i].size(); ++j) {
            result[i][j] = (indm[i][j] != 0);
        }
    }
    return result;
}

template <typename Matrix>
void print_matrix(Matrix const &m)
{
    for (auto const &row : m)
    {
        using T = typename iterator_traits<decltype(begin(row))>::value_type; //iterator_traits allowed me to get rid of the print_matrix specialization of C-style array
        copy(begin(row), end(row), ostream_iterator<T>(cout, " "));
        cout << "\n";
    }
}


int main()
{
    Graph g(11);
    g.addEdge(0, 1);
    g.addEdge(0, 6);
    g.addEdge(0, 7);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(2, 7);
    g.addEdge(3, 4);
    g.addEdge(4, 5);
    g.addEdge(4, 8);
    g.addEdge(5, 6);
    g.addEdge(7, 8);
    g.addEdge(8, 9);
    g.addEdge(9, 10);
    g.addEdge(10, 6);

    // std::vector<int> cycles;
    //g.Gotlieb(back_inserter(cycles));
    vector<vector<int>> treeAdjMat = g.getTreeAdjMat();
    // std::ofstream of("cycles.data");
    // print_cycles(begin(cycles), end(cycles), cout);
    // print_cycles(begin(cycles), end(cycles), of);

    // g.printMat();
    // print_matrix(treeAdjMat);
    // cout << "\n";

    std::istringstream stream(circuit_matrix_text);
    vector<vector<int>> matrix;

    std::string line;
    while (std::getline(stream, line)) {
        std::istringstream line_stream(line);

        int element;
        vector<int> row;

        while (line_stream >> element) {
            row.push_back(element);
        }

        matrix.push_back(row);
    }

    vector<vector<bool>> matrixBool = makeBoolMatrixFromIntMatrix(matrix);
    print_matrix(matrixBool);
    int r = matrixBool.size();
    int m;
    vector<vector<int>> output = g.Gotlieb(r, &m, matrixBool);
    printIndm(output);
    std::ofstream outfile1("indm.txt");
    printIndmTo(outfile1, output);
}

