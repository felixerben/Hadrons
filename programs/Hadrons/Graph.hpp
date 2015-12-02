/*
 * Graph.hpp, part of Grid
 *
 * Copyright (C) 2015 Antonin Portelli
 *
 * Grid is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Grid is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Grid.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef Hadrons_Graph_hpp_
#define Hadrons_Graph_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                          Oriented graph class                              *
 ******************************************************************************/
// I/O for edges
template <typename T>
std::ostream & operator<<(std::ostream &out, const std::pair<T, T> &e)
{
    out << "\""  << e.first << "\" -> \"" << e.second << "\"";
    
    return out;
}

// main class
template <typename T>
class Graph
{
public:
    typedef std::pair<T, T> Edge;
public:
    // constructor
    Graph(void) = default;
    // destructor
    virtual ~Graph(void) = default;
    // access
    void         addVertex(const T &value);
    void         addEdge(const Edge &e);
    void         addEdge(const T &start, const T &end);
    void         removeVertex(const T &value);
    void         removeEdge(const Edge &e);
    void         removeEdge(const T &start, const T &end);
    unsigned int size(void) const;
    // tests
    bool gotValue(const T &value) const;
    // graph topological manipulations
    std::vector<T>        getAdjacentVertices(const T &value) const;
    std::vector<T>        getChildren(const T &value) const;
    std::vector<T>        getParents(const T &value) const;
    std::vector<T>        getRoots(void) const;
    std::vector<Graph<T>> getConnectedComponents(void) const;
    std::stack<T>         topoSort(void);
    // I/O
    friend std::ostream & operator<<(std::ostream &out, const Graph<T> &g)
    {
        out << "{";
        for (auto &e: g.edgeSet_)
        {
            out << e << ", ";
        }
        if (g.edgeSet_.size() != 0)
        {
            out << "\b\b";
        }
        out << "}";
        
        return out;
    }
private:
    // vertex marking
    void      mark(const T &value, const bool doMark = true);
    void      markAll(const bool doMark = true);
    void      unmark(const T &value);
    void      unmarkAll(void);
    bool      isMarked(const T &value) const;
    const T * getFirstMarked(const bool isMarked = true) const;
    const T * getFirstUnmarked(void) const;
    // prune marked/unmarked vertices
    void removeMarked(const bool isMarked = true);
    void removeUnmarked(void);
    // depth-first search marking
    void depthFirstSearch(void);
    void depthFirstSearch(const T &root);
private:
    std::map<T, bool> isMarked_;
    std::set<Edge>    edgeSet_;
};

/******************************************************************************
 *                       template implementation                              *
 ******************************************************************************/
// access //////////////////////////////////////////////////////////////////////
template <typename T>
void Graph<T>::addVertex(const T &value)
{
    isMarked_[value] = false;
}

template <typename T>
void Graph<T>::addEdge(const Edge &e)
{
    addVertex(e.first);
    addVertex(e.second);
    edgeSet_.insert(e);
}

template <typename T>
void Graph<T>::addEdge(const T &start, const T &end)
{
    addEdge(Edge(start, end));
}

template <typename T>
void Graph<T>::removeVertex(const T &value)
{
    // remove vertex from the mark table
    auto vIt = isMarked_.find(value);
    
    if (vIt != isMarked_.end())
    {
        isMarked_.erase(vIt);
    }
    else
    {
        HADRON_ERROR("vertex " << value << " does not exists");
    }

    // remove all edges containing the vertex
    auto pred = [&value](const Edge &e)
    {
        return ((e.first == value) or (e.second == value));
    };
    auto eIt = find_if(edgeSet_.begin(), edgeSet_.end(), pred);
    
    while (eIt != edgeSet_.end())
    {
        edgeSet_.erase(eIt);
        eIt = find_if(edgeSet_.begin(), edgeSet_.end(), pred);
    }
}

template <typename T>
void Graph<T>::removeEdge(const Edge &e)
{
    auto eIt = edgeSet_.find(e);
    
    if (eIt != edgeSet_.end())
    {
        edgeSet_.erase(eIt);
    }
    else
    {
        HADRON_ERROR("edge "  << e << " does not exists");
    }
}

template <typename T>
void Graph<T>::removeEdge(const T &start, const T &end)
{
    removeEdge(Edge(start, end));
}

template <typename T>
unsigned int Graph<T>::size(void) const
{
    return isMarked_.size();
}

// tests ///////////////////////////////////////////////////////////////////////
template <typename T>
bool Graph<T>::gotValue(const T &value) const
{
    try
    {
        isMarked_.at(value);
    }
    catch (std::out_of_range &)
    {
        return false;
    }

    return true;
}

// vertex marking //////////////////////////////////////////////////////////////
template <typename T>
void Graph<T>::mark(const T &value, const bool doMark)
{
    try
    {
        isMarked_.at(value) = doMark;
    }
    catch (std::out_of_range &)
    {
        HADRON_ERROR("vertex " << value << " does not exists");
    }
}

template <typename T>
void Graph<T>::markAll(const bool doMark)
{
    for (auto &v: isMarked_)
    {
        mark(v.first, doMark);
    }
}

template <typename T>
void Graph<T>::unmark(const T &value)
{
    mark(value, false);
}

template <typename T>
void Graph<T>::unmarkAll(void)
{
    markAll(false);
}

template <typename T>
bool Graph<T>::isMarked(const T &value) const
{
    try
    {
        return isMarked_.at(value);
    }
    catch (std::out_of_range &)
    {
        HADRON_ERROR("vertex " << value << " does not exists");
        
        return false;
    }
}

template <typename T>
const T * Graph<T>::getFirstMarked(const bool isMarked) const
{
    auto pred = [&isMarked](const std::pair<T, bool> &v)
    {
        return (v.second == isMarked);
    };
    auto vIt = std::find_if(isMarked_.begin(), isMarked_.end(), pred);
    
    if (vIt != isMarked_.end())
    {
        return &(vIt->first);
    }
    else
    {
        return nullptr;
    }
}

template <typename T>
const T * Graph<T>::getFirstUnmarked(void) const
{
    return getFirstMarked(false);
}

// prune marked/unmarked vertices //////////////////////////////////////////////
template <typename T>
void Graph<T>::removeMarked(const bool isMarked)
{
    auto isMarkedCopy = isMarked_;
    
    for (auto &v: isMarkedCopy)
    {
        if (v.second == isMarked)
        {
            removeVertex(v.first);
        }
    }
}

template <typename T>
void Graph<T>::removeUnmarked(void)
{
    removeMarked(false);
}

// depth-first search marking //////////////////////////////////////////////////
template <typename T>
void Graph<T>::depthFirstSearch(void)
{
    depthFirstSearch(isMarked_.begin()->first);
}

template <typename T>
void Graph<T>::depthFirstSearch(const T &root)
{
    std::vector<T> adjacentVertex;
    
    mark(root);
    adjacentVertex = getAdjacentVertices(root);
    for (auto &v: adjacentVertex)
    {
        if (!isMarked(v))
        {
            depthFirstSearch(v);
        }
    }
}

// graph topological manipulations /////////////////////////////////////////////
template <typename T>
std::vector<T> Graph<T>::getAdjacentVertices(const T &value) const
{
    std::vector<T> adjacentVertex;
    
    auto pred = [&value](const Edge &e)
    {
        return ((e.first == value) or (e.second == value));
    };
    auto eIt = find_if(edgeSet_.begin(), edgeSet_.end(), pred);
    
    while (eIt != edgeSet_.end())
    {
        if (eIt->first == value)
        {
            adjacentVertex.push_back((*eIt).second);
        }
        else if (eIt->second == value)
        {
            adjacentVertex.push_back((*eIt).first);
        }
        eIt = find_if(++eIt, edgeSet_.end(), pred);
    }
    
    return adjacentVertex;
}

template <typename T>
std::vector<T> Graph<T>::getChildren(const T &value) const
{
    std::vector<T> child;
    
    auto pred = [&value](const Edge &e)
    {
        return (e.first == value);
    };
    auto eIt = find_if(edgeSet_.begin(), edgeSet_.end(), pred);
    
    while (eIt != edgeSet_.end())
    {
        child.push_back((*eIt).second);
        eIt = find_if(++eIt, edgeSet_.end(), pred);
    }
    
    return child;
}

template <typename T>
std::vector<T> Graph<T>::getParents(const T &value) const
{
    std::vector<T> parent;
    
    auto pred = [&value](const Edge &e)
    {
        return (e.second == value);
    };
    auto eIt = find_if(edgeSet_.begin(), edgeSet_.end(), pred);
    
    while (eIt != edgeSet_.end())
    {
        parent.push_back((*eIt).first);
        eIt = find_if(++eIt, edgeSet_.end(), pred);
    }
    
    return parent;
}

template <typename T>
std::vector<T> Graph<T>::getRoots(void) const
{
    std::vector<T> root;
    
    for (auto &v: isMarked_)
    {
        auto parent = getParents(v.first);
        
        if (parent.size() == 0)
        {
            root.push_back(v.first);
        }
    }
    
    return root;
}

template <typename T>
std::vector<Graph<T>> Graph<T>::getConnectedComponents(void) const
{
    std::vector<Graph<T>> res;
    Graph<T>              copy(*this);
    
    while (copy.size() > 0)
    {
        copy.depthFirstSearch();
        res.push_back(copy);
        res.back().removeUnmarked();
        res.back().unmarkAll();
        copy.removeMarked();
        copy.unmarkAll();
    }
    
    return res;
}

// topological sort using Kahn's algorithm
template <typename T>
std::stack<T> Graph<T>::topoSort(void)
{
    std::stack<T>     res;
    const T           *vPt;
    std::map<T, bool> tmpMarked(isMarked_);
    
    // visit function
    std::function<void(const T &)> visit = [&](const T &v)
    {
        if (tmpMarked.at(v))
        {
            HADRON_ERROR("cannot topologically sort a cyclic graph");
        }
        if (!this->isMarked(v))
        {
            std::vector<T> child = this->getChildren(v);
            
            tmpMarked[v] = true;
            for (auto &c: child)
            {
                visit(c);
            }
            this->mark(v);
            tmpMarked[v] = false;
            res.push(v);
        }
    };
    
    // reset temporary marks
    for (auto &v: tmpMarked)
    {
        tmpMarked.at(v.first) = false;
    }
    // loop on unmarked vertices
    vPt = getFirstUnmarked();
    while (vPt)
    {
        visit(*vPt);
        vPt = getFirstUnmarked();
    }
    unmarkAll();
    
    return res;
}

END_HADRONS_NAMESPACE

#endif // Hadrons_Graph_hpp_