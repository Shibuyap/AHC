#include <algorithm>
#include <array>
#include <bitset>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

class SudokuSolver
{
private:
    static const int SIZE = 9;
    static const int BOX_SIZE = 3;
    array<array<int, SIZE>, SIZE> grid;
    array<bitset<SIZE + 1>, SIZE> rowUsed;
    array<bitset<SIZE + 1>, SIZE> colUsed;
    array<bitset<SIZE + 1>, SIZE> boxUsed;

    int getBoxIndex(int row, int col)
    {
        return (row / BOX_SIZE) * BOX_SIZE + (col / BOX_SIZE);
    }

    bool isValid(int row, int col, int num)
    {
        return !rowUsed[row][num] && !colUsed[col][num] && !boxUsed[getBoxIndex(row, col)][num];
    }

    void setNumber(int row, int col, int num)
    {
        grid[row][col] = num;
        rowUsed[row][num] = true;
        colUsed[col][num] = true;
        boxUsed[getBoxIndex(row, col)][num] = true;
    }

    void removeNumber(int row, int col, int num)
    {
        grid[row][col] = 0;
        rowUsed[row][num] = false;
        colUsed[col][num] = false;
        boxUsed[getBoxIndex(row, col)][num] = false;
    }

    bool solveBacktrack()
    {
        for (int row = 0; row < SIZE; row++)
        {
            for (int col = 0; col < SIZE; col++)
            {
                if (grid[row][col] == 0)
                {
                    for (int num = 1; num <= SIZE; num++)
                    {
                        if (isValid(row, col, num))
                        {
                            setNumber(row, col, num);
                            if (solveBacktrack())
                            {
                                return true;
                            }
                            removeNumber(row, col, num);

                           for(int )
                        }
                    }
                    return false;
                }
            }
        }
        return true;
    }

    bool solveConstraintPropagation()
    {
        bool changed = true;
        while (changed)
        {
            changed = false;

            // Hidden singles: 行、�?�、�?�ックス�?でそ�?�数字が入る�?�所�?1つしかな�?
            for (int num = 1; num <= SIZE; num++)
            {
                // 行での�?れたシングル
                for (int row = 0; row < SIZE; row++)
                {
                    if (!rowUsed[row][num])
                    {
                        vector<int> candidates;
                        for (int col = 0; col < SIZE; col++)
                        {
                            if (grid[row][col] == 0 && isValid(row, col, num))
                            {
                                candidates.push_back(col);
                            }
                        }
                        if (candidates.size() == 1)
                        {
                            setNumber(row, candidates[0], num);
                            changed = true;
                        }
                    }
                }

                // 列での�?れたシングル
                for (int col = 0; col < SIZE; col++)
                {
                    if (!colUsed[col][num])
                    {
                        vector<int> candidates;
                        for (int row = 0; row < SIZE; row++)
                        {
                            if (grid[row][col] == 0 && isValid(row, col, num))
                            {
                                candidates.push_back(row);
                            }
                        }
                        if (candidates.size() == 1)
                        {
                            setNumber(candidates[0], col, num);
                            changed = true;
                        }
                    }
                }

                // ボックスでの�?れたシングル
                for (int box = 0; box < SIZE; box++)
                {
                    if (!boxUsed[box][num])
                    {
                        vector<pair<int, int>> candidates;
                        int startRow = (box / BOX_SIZE) * BOX_SIZE;
                        int startCol = (box % BOX_SIZE) * BOX_SIZE;
                        for (int r = 0; r < BOX_SIZE; r++)
                        {
                            for (int c = 0; c < BOX_SIZE; c++)
                            {
                                int row = startRow + r;
                                int col = startCol + c;
                                if (grid[row][col] == 0 && isValid(row, col, num))
                                {
                                    candidates.push_back({row, col});
                                }
                            }
                        }
                        if (candidates.size() == 1)
                        {
                            setNumber(candidates[0].first, candidates[0].second, num);
                            changed = true;
                        }
                    }
                }
            }

            // Naked singles: セルに入る可能性のある数字が1つしかな�?
            for (int row = 0; row < SIZE; row++)
            {
                for (int col = 0; col < SIZE; col++)
                {
                    if (grid[row][col] == 0)
                    {
                        vector<int> candidates;
                        for (int num = 1; num <= SIZE; num++)
                        {
                            if (isValid(row, col, num))
                            {
                                candidates.push_back(num);
                            }
                        }
                        if (candidates.size() == 1)
                        {
                            setNumber(row, col, candidates[0]);
                            changed = true;
                        }
                    }
                }
            }
        }
        return isComplete();
    }

public:
    SudokuSolver()
    {
        clear();
    }

    void clear()
    {
        for (int i = 0; i < SIZE; i++)
        {
            for (int j = 0; j < SIZE; j++)
            {
                grid[i][j] = 0;
            }
            rowUsed[i].reset();
            colUsed[i].reset();
            boxUsed[i].reset();
        }
    }

    bool loadFromFile(const string &filename)
    {
        ifstream file(filename);
        if (!file.is_open())
        {
            return false;
        }

        clear();
        for (int row = 0; row < SIZE; row++)
        {
            for (int col = 0; col < SIZE; col++)
            {
                int num;
                if (!(file >> num))
                {
                    return false;
                }
                if (num != 0)
                {
                    if (num < 1 || num > 9 || !isValid(row, col, num))
                    {
                        return false;
                    }
                    setNumber(row, col, num);
                }
            }
        }
        return true;
    }

    bool loadFromInput()
    {
        clear();
        cout << "数独パズルを�?�力してください (0は空のセル):" << endl;
        for (int row = 0; row < SIZE; row++)
        {
            for (int col = 0; col < SIZE; col++)
            {
                int num;
                if (!(cin >> num))
                {
                    return false;
                }
                if (num != 0)
                {
                    if (num < 1 || num > 9 || !isValid(row, col, num))
                    {
                        cout << "無効な入力で�?: �?" << (row + 1) << ", �?" << (col + 1) << endl;
                        return false;
                    }
                    setNumber(row, col, num);
                }
            }
        }
        return true;
    }

    bool isComplete()
    {
        for (int row = 0; row < SIZE; row++)
        {
            for (int col = 0; col < SIZE; col++)
            {
                if (grid[row][col] == 0)
                {
                    return false;
                }
            }
        }
        return true;
    }

    bool solve()
    {
        auto start = chrono::high_resolution_clock::now();

        // まず制�?伝播で可能な限り解�?
        if (solveConstraintPropagation())
        {
            auto end = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
            cout << "制�?伝播のみで解決しました (時間: " << duration.count() << " マイクロ�?)" << endl;
            return true;
        }

        // バックトラ�?キングで残りを解�?
        bool result = solveBacktrack();
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end - start);

        if (result)
        {
            cout << "バックトラ�?キングで解決しました (時間: " << duration.count() << " マイクロ�?)" << endl;
        }
        else
        {
            cout << "解が見つかりませんでした (時間: " << duration.count() << " マイクロ�?)" << endl;
        }

        return result;
    }

    void printGrid()
    {
        cout << "┌───────┬───────┬───────�?" << endl;
        for (int row = 0; row < SIZE; row++)
        {
            cout << "�? ";
            for (int col = 0; col < SIZE; col++)
            {
                if (grid[row][col] == 0)
                {
                    cout << "·";
                }
                else
                {
                    cout << grid[row][col];
                }
                cout << " ";
                if (col % 3 == 2)
                {
                    cout << "�? ";
                }
            }
            cout << endl;
            if (row % 3 == 2 && row != SIZE - 1)
            {
                cout << "├───────┼───────┼───────┤" << endl;
            }
        }
        cout << "└───────┴───────┴───────�?" << endl;
    }

    bool saveToFile(const string &filename)
    {
        ofstream file(filename);
        if (!file.is_open())
        {
            return false;
        }

        for (int row = 0; row < SIZE; row++)
        {
            for (int col = 0; col < SIZE; col++)
            {
                file << grid[row][col];
                if (col < SIZE - 1)
                    file << " ";
            }
            file << endl;
        }
        return true;
    }

    int countEmptyCells()
    {
        int count = 0;
        for (int row = 0; row < SIZE; row++)
        {
            for (int col = 0; col < SIZE; col++)
            {
                if (grid[row][col] == 0)
                {
                    count++;
                }
            }
        }
        return count;
    }

    bool isValidSudoku()
    {
        for (int row = 0; row < SIZE; row++)
        {
            for (int col = 0; col < SIZE; col++)
            {
                if (grid[row][col] != 0)
                {
                    int num = grid[row][col];
                    removeNumber(row, col, num);
                    if (!isValid(row, col, num))
                    {
                        setNumber(row, col, num);
                        return false;
                    }
                    setNumber(row, col, num);
                }
            }
        }
        return true;
    }
};

int main()
{
    SudokuSolver solver;
    int choice;

    cout << "=== 数独ソルバ�?� ===" << endl;
    cout << "1. ファイルから読み込み" << endl;
    cout << "2. 手動入�?" << endl;
    cout << "3. �?ストケース実�?" << endl;
    cout << "選択してください (1-3): ";
    cin >> choice;

    switch (choice)
    {
    case 1:
    {
        string filename;
        cout << "ファイル名を入力してください: ";
        cin >> filename;

        if (!solver.loadFromFile("./in/" + filename))
        {
            cout << "ファイルの読み込みに失敗しました�?" << endl;
            return 1;
        }
        break;
    }
    case 2:
    {
        if (!solver.loadFromInput())
        {
            cout << "入力に失敗しました�?" << endl;
            return 1;
        }
        break;
    }
    case 3:
    {
        // �?ストケースを作�??
        cout << "�?ストケースを実行しま�?..." << endl;
        if (!solver.loadFromFile("./in/test1.txt"))
        {
            cout << "�?ストファイルが見つかりません。サンプルパズルを使用します�?" << endl;
            // 中級レベルのサンプルパズル
            solver.clear();
        }
        break;
    }
    default:
        cout << "無効な選択です�?" << endl;
        return 1;
    }

    if (!solver.isValidSudoku())
    {
        cout << "無効な数独パズルです�?" << endl;
        return 1;
    }

    cout << "\n初期状�?:" << endl;
    solver.printGrid();
    cout << "空のセル数: " << solver.countEmptyCells() << endl;

    cout << "\n解析中..." << endl;
    if (solver.solve())
    {
        cout << "\n解�?:" << endl;
        solver.printGrid();

        string save;
        cout << "\n解答をファイルに保存しますか? (y/n): ";
        cin >> save;
        if (save == "y" || save == "Y")
        {
            string outfile;
            cout << "出力ファイル名を入力してください: ";
            cin >> outfile;
            if (solver.saveToFile("./out/" + outfile))
            {
                cout << "保存しました: ./out/" << outfile << endl;
            }
            else
            {
                cout << "保存に失敗しました�?" << endl;
            }
        }
    }
    else
    {
        cout << "こ�?�数独は解けません�?" << endl;
    }

    return 0;
}
