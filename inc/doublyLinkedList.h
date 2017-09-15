/** 
 * @file doublyLinkedList.h
 * @brief Implements DoublyLinkedList.
 * @author Markus Schmidt
 */
#ifndef DOUBLY_LINKED_LIST
#define DOUBLY_LINKED_LIST


template <typename Content>
class DoublyLinkedList;

/**
 * @brief The elements of the doubly linked list.
 */
template <typename Content>
class DoublyLinkedListElement
{
private:
	/* pointer to the the last element of the linked list */
	std::weak_ptr<DoublyLinkedListElement<Content>> pxLastNode;
	/* the content of the list element */
	std::shared_ptr<Content> pxContent;
	/* pointer to the next element of the linked list */
	std::shared_ptr<DoublyLinkedListElement<Content>> pxNextNode;

	//Deadlock prevention technique: always lock the elements IN-order of the list
	/* when multi threading the nodes need to be locked before inserting removing elements of the list */
	std::mutex xMutex;

public:
	/**
	 * @brief Constructs a new List Element and inserts it between the two given nodes.
	 */
	DoublyLinkedListElement(std::weak_ptr<DoublyLinkedListElement<Content>> pxLastNode,
		std::shared_ptr<Content> pxContent, std::shared_ptr<DoublyLinkedListElement<Content>> pxNextNode)
		:
		pxLastNode(pxLastNode),
		pxContent(pxContent),
		pxNextNode(pxNextNode)
	{}//constructor
	/**
	 * @brief Constructs a new List Element.
	 */
	DoublyLinkedListElement(){}//constructor

	/**
	 * @brief sets the next pointer to the given node.
	 * @note not thread save!
	 */
	void setNextNode(std::shared_ptr<DoublyLinkedListElement<Content>> pxNextNode)
	{
		this->pxNextNode = pxNextNode;
	}//function

	/**
	 * @brief sets the last pointer to the given node.
	 * @note not thread save!
	 */
	void setLastNode(std::weak_ptr<DoublyLinkedListElement<Content>> pxLastNode)
	{
		this->pxLastNode = pxLastNode;
	}//function
	
	/**
	 * @brief The next Node.
	 */
	const std::shared_ptr<DoublyLinkedListElement<Content>> getNextNode() const 
	{ 
		return pxNextNode;
	}//function
	
	/**
	 * @brief The last Node.
	 */
	const std::shared_ptr<DoublyLinkedListElement<Content>> getLastNode() const 
	{ 
		return pxLastNode.lock();
	}//function
	
	/**
	 * @brief The next Node.
	 */
	std::shared_ptr<DoublyLinkedListElement<Content>> getNextNode() 
	{ 
		return pxNextNode; 
	}//function
	
	/**
	 * @brief The last Node.
	 */
	std::shared_ptr<DoublyLinkedListElement<Content>> getLastNode() 
	{ 
		return pxLastNode.lock(); 
	}//function

	/**
	 * @brief Mutex sed to lock a list element.
	 * @note The user of the linked list will never have access to an object of 
	 * DoublyLinkedListElement so the mutex will be used from DoublyLinkedList only.
	 */
	std::mutex *getMutex() { return &xMutex; }//function

	/**
	 * @brief returns true.
	 * @details
	 * The two ends of the doubly linked list are object from the class DoublyLinkedListEnd. 
	 * This method is used to distinguish between DoublyLinkedListEnd and DoublyLinkedListElement.
	 */
	virtual bool isListElement() const { return true; }//function

	/**
	 * @brief Content of the list.
	 */
	virtual std::shared_ptr<Content> getContent() const { return pxContent; }//function
};//class

/**
 * @brief Both ends of the list consist of DoublyLinkedListEnd objects.
 */
template <typename Content>
class DoublyLinkedListEnd : public DoublyLinkedListElement<Content>
{

public:
	DoublyLinkedListEnd()
		:
		DoublyLinkedListElement<Content>()
	{}

	/**
	 * @brief returns false.
	 * @details
	 * The two ends of the doubly linked list are object from the class DoublyLinkedListEnd. 
	 * This method is used to distinguish between DoublyLinkedListEnd and DoublyLinkedListElement.
	 */
	bool isListElement() const { return false; };//function

	/**
	 * @brief returns nullptr.
	 */
	std::shared_ptr<Content> getContent() const { return nullptr; }//function
};//class

/**
 * @brief The doubly linked list.
 */
template <typename Content>
class DoublyLinkedList
{

private:
	/** points to the left leaf of the list; used to determine the beginning of the list.*/
	std::shared_ptr<DoublyLinkedListEnd<Content>> pxFrirstLeaf = nullptr;
	/** points to the right leaf of the list; used to determine the end of the list. */
	std::shared_ptr<DoublyLinkedListEnd<Content>> pxLastLeaf = nullptr;

	/** returns the first element (containing data) of the list;*/
	const std::shared_ptr<DoublyLinkedListElement<Content>> getRoot() const 
	{ 
		return pxFrirstLeaf->getNextNode(); 
	}

	/** returns the first element (containing data) of the list;*/
	std::shared_ptr<DoublyLinkedListElement<Content>> getRoot() 
	{ 
		return pxFrirstLeaf->getNextNode(); 
	}
public:

	/**
	 * @brief A list Interator.
	 * @details
	 * Objects of Iterator are used to give any users of the list access to the 
	 * individual list elements.
	*/
	class Iterator{
		/* 
		 * allow doubly linked list to access all private members
		 * necessary because we need to access the actual DoublyLinkedListElement in DoublyLinkedList
		 * but all other classes should not have access to anything but the Content
		 */
		friend class DoublyLinkedList;
	private:
		std::shared_ptr<DoublyLinkedListElement<Content>> pxPos;

		std::shared_ptr<DoublyLinkedListElement<Content>> getElement(){ return pxPos; }//function
	public:
		/**
		 * @brief Create a new iterator that points to the given list element.
		 */
		Iterator(std::shared_ptr<DoublyLinkedListElement<Content>> pxPos)
			:
			pxPos(pxPos)
		{}//constructor

		/**
		 * @brief Copy from another iterator.
		 */
		Iterator(const Iterator &rxOther)
			:
			pxPos(rxOther.pxPos)
		{}//constructor

		/**
		 * @brief Move to the next list element.
		 */
		void operator++(){ pxPos = pxPos->getNextNode(); }//function

		/**
		 * @brief Move to the previous list element.
		 */
		void operator--(){ pxPos = pxPos->getLastNode(); }//function
		
		/**
		 * @brief get a copy of this iterator.
		 */
		Iterator getCopy() const { return Iterator(*this); }//function
		
		/**
		 * @brief Copy from another iterator.
		 */
		void operator=(const Iterator &rxOther){ pxPos = rxOther.pxPos; }//function
		
		/**
		 * @brief Get the content of the iterator.
		 */
		std::shared_ptr<Content> operator*() const { return pxPos->getContent(); }//function

		/**
		 * @brief Get the content of the iterator.
		 */
		Content* operator->() const { return pxPos->getContent().get(); }//function
		
		/**
		 * @brief Returns true if the iterator is actually sitting on a list element at the moment.
		 */
		bool isListElement() const { return pxPos->isListElement(); }//function
	};

	/**
	* @brief Sets up an empty doubly linked list.
	*/
	DoublyLinkedList()
		:
		pxFrirstLeaf(new DoublyLinkedListEnd<Content>()),
		pxLastLeaf(new DoublyLinkedListEnd<Content>())
	{
		pxFrirstLeaf->setNextNode(pxLastLeaf);
		pxLastLeaf->setLastNode(pxFrirstLeaf);
	}//constructor

	/**
	* @brief Get an Iterator pointing to the beginning of the list.
	*/
	Iterator begin() const { return Iterator(getRoot()); }//function
	
	/**
	* @brief Get an Iterator pointing to the end of the list.
	*/
	Iterator end() const { return Iterator(pxLastLeaf->getLastNode()); }//function

	
	/**
	* @brief Execute the given function for each element in the list.
	*/
	void forEach(std::function<void(std::shared_ptr<Content>)> fDo) const
	{
		std::shared_ptr<DoublyLinkedListElement<Content>> pxCurr = getRoot();
		while (pxCurr->isListElement())
		{
			fDo(pxCurr->getContent());
			pxCurr = pxCurr->getNextNode();
		}//while
	}//function

	/**
	* @brief Remove the list element pointed to by the iterator.
	* @detail
	* Invalidates the iterator.
	*/
	void removeNode(Iterator pxItt)
	{
		removeNode(pxItt.getElement());
	}//function

	//boost python does not like overloaded functions.
	//create a un overloaded function for boost
	///@brief for boost python
	void removeNode_boost(Iterator pxItt)
	{
		removeNode(pxItt.getElement());
	}//function

	/**
	* @brief Remove the give list element.
	*/
	void removeNode(std::shared_ptr<DoublyLinkedListElement<Content>> pxNode)
	{
		//lock the last, this and the next node
		//Deadlock prevention technique: always lock the elements IN-order of the list
		std::lock_guard<std::mutex> xOtherGuard(*pxNode->getLastNode()->getMutex());
		std::lock_guard<std::mutex> xMyGuard(*pxNode->getMutex());
		std::lock_guard<std::mutex> xOtherGuard2(*pxNode->getNextNode()->getMutex());
		//check if this node is the root node; if so move the root node to the next node
		pxNode->getLastNode()->setNextNode(pxNode->getNextNode());
		pxNode->getNextNode()->setLastNode(pxNode->getLastNode());
	}//function

	/**
	* @brief Insert before the given list element.
	* @returns a iterator pointing to the new element.
	*/
	Iterator insertBefore(std::shared_ptr<Content> pxContent, std::shared_ptr<DoublyLinkedListElement<Content>> pxNode)
	{
		//lock the last and this node
		//Deadlock prevention technique: always lock the elements IN-order of the list
		std::lock_guard<std::mutex> xOtherGuard(*pxNode->getLastNode()->getMutex());
		std::lock_guard<std::mutex> xMyGuard(*pxNode->getMutex());

		std::shared_ptr<DoublyLinkedListElement<Content>> pxNewElement
			= std::shared_ptr<DoublyLinkedListElement<Content>>(new DoublyLinkedListElement<Content>(pxNode->getLastNode(), pxContent, pxNode));

		pxNode->getLastNode()->setNextNode(pxNewElement);
		pxNode->setLastNode(pxNewElement);

		return Iterator(pxNewElement);
	}//function

	/**
	* @brief Insert before the list element pointed to by the iterator.
	* @returns a iterator pointing to the new element.
	*/
	Iterator insertBefore(std::shared_ptr<Content> pxContent, Iterator pxItt)
	{
		return insertBefore(pxContent, pxItt.getElement());
	}//function

	//boost python does not like overloaded functions.
	//create a un overloaded function for boost
	///@brief for boost python
	Iterator insertBefore_boost(std::shared_ptr<Content> pxContent, Iterator pxItt)
	{
		return insertBefore(pxContent, pxItt.getElement());
	}//function

	/**
	* @brief Insert after the given list element.
	* @returns a iterator pointing to the new element.
	*/
	Iterator insertAfter(std::shared_ptr<Content> pxContent, std::shared_ptr<DoublyLinkedListElement<Content>> pxNode)
	{
		//lock the next and this node
		//Deadlock prevention technique: always lock the elements IN-order of the list
		std::lock_guard<std::mutex> xMyGuard(*pxNode->getMutex());
		std::lock_guard<std::mutex> xOtherGuard2(*pxNode->getNextNode()->getMutex());

		std::shared_ptr<DoublyLinkedListElement<Content>> pxNewElement
			= std::shared_ptr<DoublyLinkedListElement<Content>>(new DoublyLinkedListElement<Content>(pxNode, pxContent, pxNode->getNextNode()));

		pxNode->getNextNode()->setLastNode(pxNewElement);
		pxNode->setNextNode(pxNewElement);
			
		return Iterator(pxNewElement);
	}//function

	/**
	* @brief Insert after the list element pointed to by the iterator.
	* @returns a iterator pointing to the new element.
	*/
	Iterator insertAfter(std::shared_ptr<Content> pxContent, Iterator pxItt)
	{
		return insertAfter(pxContent, pxItt.getElement());
	}//function

	//boost python does not like overloaded functions.
	//create a un overloaded function for boost
	///@brief for boost python
	Iterator insertAfter_boost(std::shared_ptr<Content> pxContent, Iterator pxItt)
	{
		return insertAfter(pxContent, pxItt.getElement());
	}//function

	/**
	* @brief Insert at the end of the list.
	* @returns a iterator pointing to the new element.
	*/
	Iterator push_back(std::shared_ptr<Content> pxContent)
	{
		return insertBefore(pxContent, pxLastLeaf);
	}//function
	
	/**
	* @brief Insert at the front of the list.
	* @returns a iterator pointing to the new element.
	*/
	Iterator push_front(std::shared_ptr<Content> pxContent)
	{
		return insertAfter(pxContent, pxFrirstLeaf);
	}//function

	/**
	 * @brief The number of elements in the list.
	 * @note Not thread save.
	 */
	unsigned int length() const
	{
		unsigned int uiRet = 0;

		forEach([&uiRet](std::shared_ptr<Content>){ uiRet++; });//forEach

		return uiRet;
	}//function

};//class

#endif