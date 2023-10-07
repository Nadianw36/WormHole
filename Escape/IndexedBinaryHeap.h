#ifndef ESCAPE_INDEXEDBINARYHEAP_H_
#define ESCAPE_INDEXEDBINARYHEAP_H_

#include <iostream>
#include <unordered_map>
#include <vector>
#include "Escape/Graph.h"

using namespace std;
namespace Escape {
    class IndexedBinaryHeap {
    public:
        // Maps keys to their positions in the binary heap
        unordered_map<VertexIdx, int> keyToPosition;
        CGraph *cg;
        
        IndexedBinaryHeap(CGraph *g) {
            // Add a dummy element to the beginning of the binary heap to simplify indexing
            heap.push_back(0);
            cg = g;
        }

        // Get the maximum value in the indexed binary heap
       VertexIdx popMax() {
        int position = 1;
        VertexIdx maxKey = positionToKey[position];
        VertexIdx lastValue = heap.back(); //removing last element
        int lastPosition = heap.size() - 1;
        heap.pop_back();
        if (!heap.empty()) {
            heap[position] = lastValue; //replace last value as max
            keyToPosition.erase(maxKey); //remove max 
            keyToPosition[positionToKey[lastPosition]] = position;
            positionToKey[position] = positionToKey[lastPosition];
            positionToKey.erase(lastPosition);
            heapifyDown(position);
        }
        return maxKey;
    }

        // Increment the value associated with a key by 1
        void increment(int key) {
            if(keyToPosition.count(key)==0){
                heap.push_back(1);
                int position = heap.size() - 1;
                keyToPosition[key] = position;
                positionToKey[position] = key;
            } else{
                int position = keyToPosition[key];
                heap[position]++;
                heapifyDown(position);
                heapifyUp(position);
            }
        }

    private:

        vector<VertexIdx> heap; 
        // Maps positions in the binary heap to their keys
        unordered_map<int, VertexIdx> positionToKey;
        // choose vertex with highest L0 degree else choose higher overall degree
        bool greaterThan(int p1, int p2){
            return heap[p1] > heap[p2];
            // if (heap[p1] != heap[p2]) return heap[p1] > heap[p2];
            // VertexIdx dp1 = cg->degree(positionToKey[p1]);
            // VertexIdx dp2 = cg->degree(positionToKey[p2]);
            // return dp1 > dp2;
        }
        // Helper function to heapify up a position in the binary heap
        void heapifyUp(int position) {
            while (position > 1 && greaterThan(position, position/2)) {
                swap(heap[position], heap[position/2]);
                updatePositions(position, position/2);
                position /= 2;
            }
        }
    
        // Helper function to heapify down a position in the binary heap
        void heapifyDown(int position) {
            while (2*position < heap.size()) {
                int child = 2*position;
                if (child+1 < heap.size() && greaterThan(child+1, child)) {
                    child++;
                }
                if (greaterThan(child, position)) {
                    swap(heap[position], heap[child]);
                    updatePositions(position, child);
                    position = child;
                } else {
                    break;
                }
            }
        }

        // Helper function to update the key-to-position and position-to-key mappings after swapping elements in the binary heap
        void updatePositions(int position1, int position2) {
            VertexIdx key1 = positionToKey[position1];
            VertexIdx key2 = positionToKey[position2];
            keyToPosition[key1] = position2;
            keyToPosition[key2] = position1;
            positionToKey[position1] = key2;
            positionToKey[position2] = key1;
        }
    };
}
 #endif
